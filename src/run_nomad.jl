"""

	nomad(eval::Function,param::nomadParameters)

-> Run NOMAD with settings defined by param and an
optimization problem defined by eval(x).

-> Display stats from NOMAD in the REPL.

-> return an object of type *nomadResults* that contains
info about the run.

# **Arguments** :

- `eval::Function`

a function of the form :

	(success,count_eval,bb_outputs)=eval(x::Vector{Number})

`bb_outputs` being a *vector{Number}* containing
the values of objective function and constraints
for a given input vector `x`. NOMAD will seak for
minimizing the objective function and keeping
constraints inferior to 0.

`success` is a *Bool* set to `false` if the evaluation failed.

`count_eval` is a *Bool* equal to `true` if the black
box evaluation counting has to be incremented. Note
that statistic sums and averages are updated only if
`count_eval` is equal to `true`.

- `param::nomadParameters`

An object of type *Parameters* of which the
attributes are the settings of the optimization
process (dimension, output types, display options,
bounds, etc.).

# **Example** :

	using NOMAD

	function eval(x)
	    f=x[1]^2+x[2]^2
	    c=1-x[1]
		success=true
	    count_eval=true
		bb_outputs = [f,c]
	    return (success,count_eval,bb_outputs)
	end

	param = nomadParameters([5,5],["OBJ","EB"])
	#=first element of bb_outputs is the objective function, second is a
		constraint treated with the Extreme Barrier method. Initial point
		of optimization is [5,5]=#

	result = nomad(eval,param)

"""
function nomad(eval::Function,param::nomadParameters;surrogate=nothing,extended_poll=nothing)

	#=
	This function first wraps eval with a julia function eval_wrap
	that takes a C-double[] as argument and returns a C-double[].
	Then it converts all param attributes into C++ variables and
	calls the C++ function cpp_runner previously defined by
	init.
	=#

	has_sgte = !isnothing(surrogate)
	has_extpoll = !isnothing(extended_poll)

	check_eval_param(eval,param,surrogate) #check consistency of nomadParameters with problem

	m=length(param.output_types)::Int64
	n=param.dimension::Int64

	#C++ wrapper for eval(x) and surrogate
	function eval_wrap(x::Ptr{Float64})

		j_x = convert_cdoublearray_to_jlvector(x,n+1)::Vector{Float64}

		if has_sgte && j_x[n+1]==1  #last coordinate of input decides if we call the surrogate or not
			(success,count_eval,bb_outputs)=surrogate(j_x[1:n]);
		else
			(success,count_eval,bb_outputs)=eval(j_x[1:n]);
		end;
		bb_outputs=convert(Vector{Float64},bb_outputs);

		return icxx"""
	    double * c_output = new double[$m+2];
	    $:(
			#converting from Vector{Float64} to C-double[]
			for j=1:m
			    icxx"c_output[$j-1]=$(bb_outputs[j]);";
			end;
			nothing
		);

		//last coordinates of c_ouput correspond to success and count_eval
		c_output[$m]=0.0;
		c_output[$m+1]=0.0;
		if ($success) {c_output[$m]=1.0;}
		if ($count_eval) {c_output[$m+1]=1.0;}

	    return c_output;
	    """
	end

	#struct containing void pointer toward eval_wrap
	evalwrap_void_ptr_struct = @cfunction($eval_wrap, Ptr{Cdouble}, (Ptr{Cdouble},))::Base.CFunction
	#void pointer toward eval_wrap
	evalwrap_void_ptr = evalwrap_void_ptr_struct.ptr::Ptr{Nothing}

	if has_extpoll

		max_signature_index=length(param.signatures)
		sign0=nomadSignature(param.input_types)
		sign0.lower_bound=param.lower_bound
		sign0.upper_bound=param.upper_bound
		sign0.granularity=param.granularity
		pushfirst!(param.signatures,sign0)

		#C++ wrapper for extpoll(x)
		function extpoll_wrap(x::Ptr{Float64})

			j_x = convert_cdoublearray_to_jlvector(x,n)::Vector{Float64};

			(extpoll_points,signature_ind)=extended_poll(j_x);
			num_pp = length(extpoll_points)
			sizes_pp = length.(extpoll_points)

			return icxx"""
		    double * c_output = new double[1+$(sum(sizes_pp)+num_pp)];
			c_output[0] = $num_pp;
		    $:(
				index=1;
				for i=1:num_pp
					s_i=signature_ind[i]
					0<=s_i<max_signature_index || error("Extended poll error : a signature index returned by extpoll(x) is wrong")
					icxx"c_output[$index]=$(signature_ind[i]);";
					for j=1:sizes_pp[i]
				    	icxx"c_output[$(index+j)]=$(extpoll_points[i][j]);";
					end;
					index+=sizes_pp[i]+1;
				end;
				nothing
			);

		    return c_output;
		    """
		end

		#struct containing void pointer toward extpoll_wrap
		extpollwrap_void_ptr_struct = @cfunction($extpoll_wrap, Ptr{Cdouble}, (Ptr{Cdouble},))::Base.CFunction
		#void pointer toward extpoll_wrap
		extpollwrap_void_ptr = extpollwrap_void_ptr_struct.ptr::Ptr{Nothing}

	else

		extpollwrap_void_ptr=icxx"void * ptr; return ptr;"

	end



	#converting param attributes into C++ variables
	c_input_types=convert_input_types(param.input_types,n)
	c_output_types=convert_output_types(param.output_types,m)
	c_display_stats=convert_string(param.display_stats)::Cstring
	c_x0=convert_x0_to_nomadpoints_list(param.x0)
	c_lower_bound=convert_vector_to_nomadpoint(param.lower_bound)::CnomadPoint
	c_upper_bound=convert_vector_to_nomadpoint(param.upper_bound)::CnomadPoint
	c_granularity=convert_vector_to_nomadpoint(param.granularity)::CnomadPoint
	c_signatures=convert_signatures(param.signatures)

	#calling cpp_runner
	c_result = @cxx cpp_runner(param.dimension,
										length(param.output_types),
										evalwrap_void_ptr,
										extpollwrap_void_ptr,
										c_input_types,
										c_output_types,
										param.display_all_eval,
										c_display_stats,
										c_x0,
										c_lower_bound,
										c_upper_bound,
										param.max_bb_eval,
										param.max_time,
										param.display_degree,
										param.LH_init,
										param.LH_iter,
										param.sgte_cost,
										c_granularity,
										param.stop_if_feasible,
										param.VNS_search,
										(param.stat_sum_target==Inf ? 0 : param.stat_sum_target),
										param.seed,
										("STAT_AVG" in param.output_types),
										("STAT_SUM" in param.output_types),
										has_sgte,
										has_extpoll,
										c_signatures,
										param.poll_trigger,
										param.relative_trigger)

	#creating nomadResults object to return
	jl_result = nomadResults(c_result,param)

	return 	jl_result

end #nomad


######################################################
		   		#CONVERSION METHODS#
######################################################

function convert_cdoublearray_to_jlvector(c_vector,size)
	jl_vector = Vector{Float64}(undef,size)
	for i=1:size
		jl_vector[i]=icxx"return *($c_vector+$i-1);"
		if abs(jl_vector[i]-round(jl_vector[i]))<eps(jl_vector[i])^0.85
			jl_vector[i]=round(jl_vector[i])
		end
	end
	return jl_vector
end

function convert_string(jl_string)
	return pointer(jl_string)
end

function convert_vector_to_nomadpoint(jl_vector)
	size = length(jl_vector)
	return icxx"""NOMAD::Double d;
				NOMAD::Point nomadpoint($size,d);
				$:(
					for i=1:size
						icxx"nomadpoint[int($i-1)]=$(jl_vector[i]);";
					end;
					nothing
				);
				return nomadpoint;"""
end

function convert_x0_to_nomadpoints_list(jl_x0)
	return icxx"""std::vector<NOMAD::Point> c_x0;
				$:(
					for x0 in jl_x0
						xZero=convert_vector_to_nomadpoint(x0);
						icxx"c_x0.push_back($xZero);";
					end;
					nothing
				);
				return c_x0;"""
end

function convert_input_types(it,n)
	return icxx"""vector<NOMAD::bb_input_type> bbit ($n);
					$:(
						for i=1:n
							if it[i]=="B"
								icxx"bbit[$i-1]=NOMAD::BINARY;";
							elseif it[i]=="I"
								icxx"bbit[$i-1]=NOMAD::INTEGER;";
							elseif it[i]=="C"
								icxx"bbit[$i-1]=NOMAD::CATEGORICAL;";
							else
								icxx"bbit[$i-1]=NOMAD::CONTINUOUS;";
							end;
						end;
							nothing
					);
					return bbit;"""
end

function convert_output_types(ot,m)
	icxx"""vector<NOMAD::bb_output_type> bbot ($m);
			$:(
				for j=1:m
					if ot[j]=="OBJ"
						icxx"bbot[$j-1]=NOMAD::OBJ;";
					elseif ot[j]=="EB"
						icxx"bbot[$j-1]=NOMAD::EB;";
					elseif ot[j] in ["PB","CSTR"]
						icxx"bbot[$j-1]=NOMAD::PB;";
					elseif ot[j] in ["PEB","PEB_P"]
						icxx"bbot[$j-1]=NOMAD::PEB_P;";
					elseif ot[j]=="PEB_E"
						icxx"bbot[$j-1]=NOMAD::PEB_E;";
					elseif ot[j] in ["F","FILTER"]
						icxx"bbot[$j-1]=NOMAD::FILTER;";
					elseif ot[j]=="CNT_EVAL"
						icxx"bbot[$j-1]=NOMAD::CNT_EVAL;";
					elseif ot[j]=="STAT_AVG"
						icxx"bbot[$j-1]=NOMAD::STAT_AVG;";
					elseif ot[j]=="STAT_SUM"
						icxx"bbot[$j-1]=NOMAD::STAT_SUM;";
					else
						icxx"bbot[$j-1]=NOMAD::UNDEFINED_BBO;";
					end;
				end;
			);
			return bbot;"""
end


function convert_signatures(sign)
	return icxx"""std::vector<NOMAD::Signature *> c_sign;
	$:(
		for s in sign
			c_s=convert_signature(s)
			icxx"c_sign.push_back($c_s);";
		end;
		nothing
	);
	return c_sign;"""
end

function convert_signature(s)

	c_input_types=convert_input_types(s.input_types,length(s.input_types))
	c_lower_bound=convert_vector_to_nomadpoint(s.lower_bound)
	c_upper_bound=convert_vector_to_nomadpoint(s.upper_bound)
	c_granularity=convert_vector_to_nomadpoint(s.granularity)::CnomadPoint

	return icxx"""NOMAD::Display out ( std::cout );
				  out.precision ( NOMAD::DISPLAY_PRECISION_STD );
				  NOMAD::Parameters p ( out );
				  p.set_DIMENSION ($(s.dimension));
				  std::cout<<"HERE0"<<std::endl;
				  NOMAD::Double mesh_update_basis = p.get_mesh_update_basis();
				  std::cout<<"HERE1/2"<<std::endl;
				  NOMAD::Double poll_update_basis = p.get_poll_update_basis();
				  std::cout<<"HERE1"<<std::endl;
				  int mesh_coarsening_exponent = p.get_mesh_coarsening_exponent();
				  std::cout<<"HERE2"<<std::endl;
				  int mesh_refining_exponent = p.get_mesh_refining_exponent();
				  std::cout<<"HERE3"<<std::endl;
				  std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> variable_groups = p.get_variable_groups();
				  NOMAD::Signature * c_s = new NOMAD::Signature ( 	$(s.dimension),
													   		$c_input_types,
													   		$c_lower_bound,
													   		$c_upper_bound,
													   		p.get_mesh_type(),
													   		p.get_anisotropic_mesh(),
													   		p.get_anisotropy_factor(),
													   		$c_granularity,
													   		p.get_initial_poll_size(),
													   		p.get_min_poll_size(),
													   		p.get_min_mesh_size(),
													   		mesh_update_basis,
													   		poll_update_basis,
													   		mesh_coarsening_exponent,
													   		mesh_refining_exponent,
													   		p.get_initial_mesh_index(),
													   		p.get_scaling(),//
													   		p.get_fixed_variables(),//
													   		p.get_periodic_variables(),
													   		variable_groups,
													   		p.out()
													   	);

					return c_s;"""
end
