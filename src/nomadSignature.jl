mutable struct nomadSignature

    dimension::Int64
    input_types::Vector{String}
    lower_bound::Vector{Float64}
    upper_bound::Vector{Float64}

    function nomadSignature(input_types::Vector{String})
        dimension=length(input_types)
        lower_bound=[]
        upper_bound=[]
        new(dimension,input_types,lower_bound,upper_bound)
    end

    function nomadSignature(s::nomadSignature)
        new(s.dimension,copy(s.input_types),copy(s.lower_bound),copy(s.upper_bound))
    end

end
