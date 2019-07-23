## Categorical variables

The following was extracted from the NOMAD user guide.

When variables can be represented by integers, but the numbers do not mean anything and
they cannot be logically ordered without further analyses, the variables can be represented
in NOMAD using categorical variables. In particular, when your problem has a number of design
variables that can vary by selecting a parameter, this parameter can be set as a categorical variable.

# Algorithm

At the end of an iteration where categorical variables are kept fixed, if no improvement has been
made, a special step occurs, the extended poll. The extended poll first calls a user-provided
procedure defining the neighborhood of categorical variables. The procedure returns a list of points
that are neighbors of the current best point (incumbent) such that categorical variables
are changed and the other variables may or may not be changed. These points are called the
extended poll points and their dimension may be different than the current best point, for
example when a categorical variable indicates the number of continuous variables.
The functions defining the problem are then evaluated at each of the extended poll points and
the objective values are compared to the current best value. If the difference between the
objective value at the current iterate and at an extended poll point is less than a parameter
called the extended poll trigger, this extended poll point is called an extended poll center and
a new MADS run is performed from this point. This run is called an extend poll descent and
occurs on meshes that cannot be reduced more than the mesh size at the beginning of the
extended poll.
If a surrogate is available, it can be used to evaluate the neighbors during the extended poll
descent. The true functions will then be evaluated only on the most promising points. With
surrogates, the extended poll costs at most the same number of true evaluations as the number
of neighbors determined by the user-provided procedure.

# Categorical variables in NOMAD.jl

In NOMAD.jl, a categorical variable is identified by setting an input type to the
value ‘C’ in nomadParameters. In addition, solving problems with categorical variables requires to define the neighbors
of the current best point. Programming the method to define the categorical variables neighborhoods relies on a function
extended_poll(x) provided as optional argument to the method nomad() ; This function takes as argument a point (the current
iterate) and registers a list of extended poll points.
In addition, each point in the algorithm possesses a signature, indicating the characteristics related to the variables
(dimension, types and bounds). The characteristics entered in the *nomadParameters* object `param` provided to nomad()
already define a signature. But if the user wants to provide points with different signatures during the extended
poll, they must be input as a vector of *nomadSignatures* attributed to param. Moreover, for each extended poll point created by extended_poll(x), the index of the signature must be indicated (index 0 for the signature implicitly defined in nomadParameters).
More precisely, extended_poll(x) must be of the form :

    (extpoll_points_list,signatures_list) = extended_poll(x::Vector{Number})

with `extpoll_points_list` a *Vector{Vector{Number}}* containing all extended poll points, and `signatures_list` a
*Vector{Int}* containing the indices of signatures corresponding to each point (defined with the same order as in
`param.signatures`). Then, NOMAD can be run with the command :

    nomad(eval::Function,param::nomadParameters;extended_poll=extended_poll::Function)

Note that a surrogate may also be provided as another optional argument.

Although the dimension of the problem may change during optimization, the starting points must
all have the same signature.
Changing the values of the categorical variables is done exclusively during the extended poll phase
by providing the neighbors of the current best point. The logic for providing the neighbors is
entirely left to the user. For this reason, it is not necessary to provide bounds for the categorical
variables whether in the initial description of the problem or when providing extended poll
point signatures.
The main parameter for mixed variable optimization is the extended poll trigger (attribute poll_trigger
of nomadParameters). It may be given as a relative value (attribute relative_trigger of nomadParameters).
The extended poll trigger is used to compare the objective values at an extended poll point y and
at the current best point x k . If f (y) < f (x k )+trigger, then y becomes an extended poll center
from which a MADS run is performed.
Finally, please note that granularity is not available if categorical variables are used.
