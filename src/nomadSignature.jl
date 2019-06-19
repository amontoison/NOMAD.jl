mutable struct nomadSignature

    dimension::Int64
    input_types::Vector{String}
    lower_bound::Vector{Float64}
    upper_bound::Vector{Float64}
    granularity::Vector{Float64}

    function nomadSignature(n)
        dimension=n
        input_types=[]
        lower_bound=[]
        upper_bound=[]
        granularity=zeros(Float64,dimension)
        new(dimension,input_types,lower_bound,upper_bound,granularity)
    end

    function nomadSignature(s::nomadSignature)
        new(s.dimension,copy(s.input_types),copy(s.lower_bound),copy(s.upper_bound),copy(s.granularity)
    end

end
