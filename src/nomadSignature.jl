"""

    nomadSignature

mutable struct defining an additional signature for extended poll points
(see Categorical variables section in documentation for more details).

# **Constructors** :

- Classic constructor :

    `s1 = nomadSignature(input_types::Vector{String})`

- Copy constructor (deepcopy):

    `s1 = nomadSignature(s2)`

# **Attributes** :

- `input_types::Vector{String}` :

A vector containing *String* objects that define the
types and dimension of this signature (the order is important) :

String  | Input type |
:-------|:-----------|
`"R"`   | Real/Continuous |
`"B"`   | Binary |
`"I"`   | Integer |
`"C"`   | Categorical |

No default value, needs to be set.

- `lower_bound::Vector{Number}` :
Lower bound for each coordinate of the state.
If empty, no bound is taken into account.
`[]` by default.

- `upper_bound::Vector{Number}` :
Upper bound for each coordinate of the state.
If empty, no bound is taken into account.
`[]` by default.

"""
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
