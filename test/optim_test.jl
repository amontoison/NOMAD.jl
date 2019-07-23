function cost1(x)
	y=(x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2
	return y
end

function cost2(x)
	y=(x[1]+2.12*x[2]-7.02)^2+(2.04*x[1]+x[2]-5.11)^2
	return y
end

function cost3(x)
	y=(x[1]^2+x[2]^2+x[3]^2)
	return y
end

function constraint(x)
	y=2-x[1]
	return y
end

function eval1(x)
    f=cost1(x)
	success=true
    count_eval=true
	bb_outputs = [f]
    return (success,count_eval,bb_outputs)
end

function eval2(x)
    f=cost1(x)
    c=constraint(x)
	success=true
    count_eval=true
	bb_outputs = [f,c]
    return (success,count_eval,bb_outputs)
end

function eval3(x)
	f=cost2(x)
	success=true
	count_eval=true
	bb_outputs = [f]
	return (success,count_eval,bb_outputs)
end

function eval4(x)
	f=cost3(x)
	c=constraint(x)
	success=true
	count_eval=true
	bb_outputs = [f,c,f,f]
	return (success,count_eval,bb_outputs)
end

function eval5(x)
	f=cost3(x)
	success=true
	count_eval=true
	bb_outputs = [f]
	return (success,count_eval,bb_outputs)
end

function eval6(x)
	x[2:end].-=1
	if x[1]==1
		f=0.5+sum(x[2:end].^2)
	elseif x[1]==0
		f=sum(x[2:end].^2)
	else
		error("x[1]=$(x[1]) is not a boolean")
	end
	success=true
	count_eval=true
	bb_outputs = [f]
	return (success,count_eval,bb_outputs)
end

function extpoll(x)
	if !(x[1] in [0,1])
		println(x[1])
		error("x[1] is not a boolean")
	end
	if length(x[2:end])==3 && x[1]==1 #if dimension=3 and objective function is first one
		p1=[1,x[3]+2,x[4]+2]
		s1_i=1 #dimension is changed to 2
		p2=[0,x[2],x[3],x[4]]
		s2_i=0 #objective function is switched
		return ([p1,p2],[s1_i,s2_i])
	else
		p1=[0,x[2]+0.3,x[3]+0.3]
		s1_i=1
		p2=[0,x[2],x[3],x[3]]
		s2_i=0
		return ([p1,p2],[s1_i,s2_i])
	end
end

param1=nomadParameters([5,5],["OBJ"])
param1.max_bb_eval=100

param2=nomadParameters(param1)
param2.output_types=["OBJ","PB"]
param2.x0=[9,9]
param2.lower_bound=[1,1]
param2.upper_bound=[10,10]
param2.max_bb_eval=50
param2.LH_init=20

param3=nomadParameters(param1)
param3.max_time=2
param3.sgte_cost=10
param3.VNS_search=true
param3.seed=-1

param4=nomadParameters([5,5,5],["OBJ","EB","STAT_SUM","STAT_AVG"])
param4.LH_iter=3
param4.display_stats="bbe ( sol ) obj ; stat_avg ; stat_sum"
param4.stat_sum_target=50000

param5=nomadParameters([5,1,5.2],["OBJ"])
param5.input_types=["I","B","R"]
param5.granularity[3]=0.2

param6=nomadParameters([[-14,70],[1,2]],["OBJ","PB"])
param6.display_all_eval=true
param6.stop_if_feasible=true

param7=nomadParameters([1,10,10,10],["OBJ"]) #first coordinate choose the objective function
param7.max_bb_eval = 200
param7.input_types=["C","R","R","R"]
sign1=nomadSignature(["C","R","R"]) #dimension can be modified
push!(param7.signatures,sign1)

#classic run
result1 = nomad(eval1,param1)
@test result1.success
test_results_consistency(result1,param1,eval1)
@test result1.best_feasible ≈ [1.0, 3.0]
disp(result1)

#PB constraint + LH initialization + bounding box
result2 = nomad(eval2,param2)
@test result2.success
test_results_consistency(result2,param2,eval2)
disp(result2)

#surrogate + VNS search
result3 = nomad(eval3,param3;surrogate=eval1) #eval1 as a surrogate of eval3
@test result3.success
test_results_consistency(result3,param3,eval3)
disp(result3)

#EB constraint + statistic sum + statistic average + LH iterations
result4 = nomad(eval4,param4)
@test result4.success
test_results_consistency(result4,param4,eval4)
disp(result4)

#Binary and integer variables + granularity
result5 = nomad(eval5,param5)
@test result5.success
test_results_consistency(result5,param5,eval5)
disp(result5)

#stop if feasible + several initial points
result6 = nomad(eval2,param6)
@test result6.success
test_results_consistency(result6,param6,eval2)
disp(result6)

#categorical variables
result7 = nomad(eval6,param7;extended_poll=extpoll)
@test result7.success
disp(result7)
