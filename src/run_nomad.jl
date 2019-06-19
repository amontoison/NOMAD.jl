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

	if !has_sgte

		#C++ wrapper for eval(x)
		function eval_wrap(x::Ptr{Float64})

			j_x = convert_cdoublearray_to_jlvector(x,n)::Vector{Float64};
			(success,count_eval,bb_outputs)=eval(j_x);
			bb_outputs=Float64.(bb_outputs);

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

	else

		#C++ wrapper for eval(x) and surrogate
		function eval_wrap(x::Ptr{Float64})

			j_x = convert_cdoublearray_to_jlvector(x,n+1)::Vector{Float64};

			if convert(Bool,j_x[n+1]) #last coordinate of input decides if we call the surrogate or not
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

	end

	#struct containing void pointer toward eval_wrap
	evalwrap_void_ptr_struct = @cfunction($eval_wrap, Ptr{Cdouble}, (Ptr{Cdouble},))::Base.CFunction
	#void pointer toward eval_wrap
	evalwrap_void_ptr = evalwrap_void_ptr_struct.ptr::Ptr{Nothing}

	if has_extpoll

		#C++ wrapper for extpoll(x)
		function extpoll_wrap(x::Ptr{Float64})

			j_x = convert_cdoublearray_to_jlvector(x,n)::Vector{Float64};
			extpoll_points=extended_poll(j_x);
			num_pp = length(extpoll_points)
			sizes_pp = length.(extpoll_points)

			return icxx"""
		    double * c_output = new double[1+$(sum(sizes_pp)+num_pp)];
			c_output[0] = $num_pp;
		    $:(
				index=1
				for i=1:num_pp
					icxx"c_output[$index]=$(sizes_pp[i]);";
					for j=1:sizes_pp[i]
				    	icxx"c_output[$(index+j)]=$(extpoll_points[i][j]);";
					end;
					index+=sizes_pp[i]+1
				end;
				nothing
			);

		    return c_output;
		    """
		end

	else

		function extpollwrap()
			return 0;
		end

	end

	#struct containing void pointer toward extpoll_wrap
	extpollwrap_void_ptr_struct = @cfunction($extpoll_wrap, Ptr{Cdouble}, (Ptr{Cdouble},))::Base.CFunction
	#void pointer toward extpoll_wrap
	extpollwrap_void_ptr = extpollwrap_void_ptr_struct.ptr::Ptr{Nothing}

	#converting param attributes into C++ variables
	c_input_types=convert_vectorstring(param.input_types,n)
	c_output_types=convert_vectorstring(param.output_types,m)
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
	end
	return jl_vector
end

function convert_vectorstring(jl_vectorstring,size)
	return icxx"""std::vector<std::string> c_vectorstring;
					$:(
						for i=1:size
							icxx"c_vectorstring.push_back($(pointer(jl_vectorstring[i])));";
						end;
						nothing
					);
					return c_vectorstring;"""
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

function convert_input_types(param.input_types,n)
	return icxx"""vector<NOMAD::bb_input_type> bbit (n);
					for (int i = 0; i < n; ++i) {
						if (input_types_[i]=="B")
							{bbit[i]=NOMAD::BINARY;}
						else if (input_types_[i]=="I")
							{bbit[i]=NOMAD::INTEGER;}
						else if (input_types_[i]=="C")
							{bbit[i]=NOMAD::CATEGORICAL;}
						else
							{bbit[i]=NOMAD::CONTINUOUS;}
					}
					return bbit;"""
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
	return icxx"""NOMAD_Signature * c_s = new Signature ( s.dimension                         ,
														 bbit_1                     ,
														 d0_1                       ,
														 lb_1                       ,
														 ub_1                       ,
														 p.get_direction_types   () ,
														 p.get_sec_poll_dir_types() ,
														 p.get_int_poll_dir_types() ,
														 _p.out()                     )
