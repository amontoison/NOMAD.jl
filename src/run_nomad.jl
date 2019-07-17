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
function nomad(eval::Function,param_::nomadParameters;surrogate=nothing,extended_poll=nothing)

	#=
	This function first wraps eval with a julia function eval_wrap
	that takes a C-double[] as argument and returns a C-double[].
	Then it converts all param attributes into C++ variables and
	calls the C++ function cpp_runner previously defined by
	init.
	=#

	param=deepcopy(param_)

	has_sgte = !isnothing(surrogate)
	has_extpoll = !isnothing(extended_poll)

	check_eval_param(eval,param,surrogate) #check consistency of nomadParameters with problem

	m=length(param.output_types)::Int64

	#C++ wrapper for eval(x) and surrogate
	function eval_wrap(x::Ptr{Float64})

		n = Int64(icxx"int n = ($x)[0];
						return n;")

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
		pushfirst!(param.signatures,sign0)

		#C++ wrapper for extpoll(x)
		function extpoll_wrap(x::Ptr{Float64})

			n = Int64(icxx"int n = ($x)[0];
							return n;")

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
					0<=s_i<=max_signature_index || error("Extended poll error : a signature index returned by extpoll(x) is wrong")
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

	c_out = icxx"""int argc;
					char ** argv;
					NOMAD::Display out ( std::cout );
					out.precision ( NOMAD::DISPLAY_PRECISION_STD );
					NOMAD::begin ( argc , argv );
					return out;"""

	c_parameter = convert_parameter(param,m,has_sgte,has_extpoll,c_out)

	c_signatures = convert_signatures(param.signatures,c_parameter)

	#calling cpp_runner
	c_result = @cxx cpp_runner(c_parameter,
								c_out,
								length(param.output_types),
								evalwrap_void_ptr,
								extpollwrap_void_ptr,
								("STAT_AVG" in param.output_types),
								("STAT_SUM" in param.output_types),
								has_sgte,
								has_extpoll,
								c_signatures)

	#creating nomadResults object to return
	jl_result = nomadResults(c_result,param)

	return 	jl_result

end #nomad



######################################################
		   		#CONVERSION METHODS#
######################################################

function convert_parameter(param,m,has_sgte,has_extpoll,out)

	#converting param attributes into C++ variables
	c_input_types=convert_input_types(param.input_types,param.dimension)
	c_output_types=convert_output_types(param.output_types,m)
	c_display_stats=convert_string(param.display_stats)
	c_x0=convert_x0_to_nomadpoints_list(param.x0)
	c_lower_bound=convert_lower_bound(param.lower_bound,param.input_types)
	c_upper_bound=convert_upper_bound(param.upper_bound,param.input_types)
	c_granularity=convert_vector_to_nomadpoint(param.granularity)

	return icxx"""NOMAD::Parameters * p = new NOMAD::Parameters( $out );

			p->set_DIMENSION ($(param.dimension));
			p->set_BB_INPUT_TYPE ( $c_input_types );
			p->set_BB_OUTPUT_TYPE ( $c_output_types );
			p->set_DISPLAY_ALL_EVAL( $(param.display_all_eval) );
			p->set_DISPLAY_STATS( $c_display_stats );
			for (int i = 0; i < ($c_x0).size(); ++i) {p->set_X0( ($c_x0)[i] );}  // starting points
			p->set_LOWER_BOUND( $c_lower_bound );
			p->set_UPPER_BOUND( $c_upper_bound );
			if ($(param.max_bb_eval)>0) {p->set_MAX_BB_EVAL($(param.max_bb_eval));}
			if ($(param.max_time)>0) {p->set_MAX_TIME($(param.max_time));}
			p->set_DISPLAY_DEGREE($(param.display_degree));
			p->set_HAS_SGTE($has_sgte);
			if ($has_sgte) {p->set_SGTE_COST($(param.sgte_cost));}
			p->set_STATS_FILE("temp.txt","bbe | sol | bbo");
			p->set_LH_SEARCH($(param.LH_init),$(param.LH_iter));
			if ($(!has_extpoll)) {p->set_GRANULARITY($c_granularity);}
			p->set_STOP_IF_FEASIBLE($(param.stop_if_feasible));
			p->set_VNS_SEARCH($(param.VNS_search));
			if ($(param.stat_sum_target)>0) {p->set_STAT_SUM_TARGET($(param.stat_sum_target));}
			p->set_SEED($(param.seed));
			p->set_EXTENDED_POLL_ENABLED($has_extpoll);
			if ($has_extpoll) {p->set_EXTENDED_POLL_TRIGGER ( $(param.poll_trigger) , $(param.relative_trigger) );}

			p->check();
			// parameters validation

			return p;"""
end

function convert_cdoublearray_to_jlvector(c_vector,size)
	jl_vector = Vector{Float64}(undef,size)
	for i=1:size
		jl_vector[i]=icxx"return *($c_vector+$i);" #coordinate 0 is not taken into account because it indicates dimension
		if abs(jl_vector[i]-round(jl_vector[i]))<eps(jl_vector[i])^0.85
			jl_vector[i]=round(jl_vector[i])
		end
	end
	return jl_vector
end

function convert_string(jl_string)
	return pointer(jl_string)
end

function convert_lower_bound(jl_lb,it)
	dim = length(jl_lb)
	return icxx"""NOMAD::Point lb ($dim);
				$:(
					for i=1:dim
						if it[i]!="C" && jl_lb[i]>-Inf
							icxx"lb[int($i-1)]=$(jl_lb[i]);";
						end;
					end;
					nothing
				);
				return lb;"""
end

function convert_upper_bound(jl_ub,it)
	dim = length(jl_ub)
	return icxx"""NOMAD::Point ub ($dim);
				$:(
					for i=1:dim
						if it[i]!="C" && jl_ub[i]<Inf
							icxx"ub[int($i-1)]=$(jl_ub[i]);";
						end;
					end;
					nothing
				);
				return ub;"""
end

function convert_vector_to_nomadpoint(jl_vector)
	dim = length(jl_vector)
	return icxx"""NOMAD::Point nomadpoint($dim);
				$:(
					for i=1:dim
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


function convert_signatures(sign,c_param)
	return icxx"""std::vector<NOMAD::Signature *> c_sign;
	$:(
		for s in sign
			c_s=convert_signature(s,c_param)
			icxx"c_sign.push_back($c_s);";
		end;
		nothing
	);
	return c_sign;"""
end

function convert_signature(s,c_param)

	c_input_types=convert_input_types(s.input_types,s.dimension)
	c_lower_bound=convert_lower_bound(s.lower_bound,s.input_types)
	c_upper_bound=convert_upper_bound(s.upper_bound,s.input_types)
	c_init_poll_size=generate_init_poll_size(s)

	return icxx"""NOMAD::Parameters * p = $c_param;
					NOMAD::Signature * c_s = new NOMAD::Signature ( $(s.dimension),
																  $c_input_types,
																  $c_init_poll_size,
																  $c_lower_bound,
																  $c_upper_bound,
																  p->get_direction_types(),
																  p->get_sec_poll_dir_types(),
																  p->get_int_poll_dir_types(),
																  p->out()
																);
	                return c_s;"""
end

function generate_init_poll_size(s)
	return icxx"""NOMAD::Point init_poll_size ( $(s.dimension) );
				$:(
					for i=1:s.dimension
						if s.input_types[i]=="C"
					  	  icxx"init_poll_size[int($i-1)]=1.0;";
						elseif s.lower_bound[i]>-Inf && s.upper_bound[i]<Inf
						  icxx"init_poll_size[int($i-1)]=$(0.1*(s.upper_bound[i]-s.lower_bound[i]));";
						elseif s.lower_bound[i]>-Inf
						  icxx"init_poll_size[int($i-1)]=$(0.1*abs(s.lower_bound[i]));";
						elseif s.upper_bound[i]<Inf
						  icxx"init_poll_size[int($i-1)]=$(0.1*abs(s.upper_bound[i]));";
						else
						  icxx"init_poll_size[int($i-1)]=1.0;";
						end;
					end;
					nothing
				);
				return init_poll_size;"""
end
