"""

	init(".../nomad.3.9.1")

load NOMAD libraries and create C++ class and function
needed to handle NOMAD optimization process.

This function has to be called once before using runopt.
It is automatically called when importing NOMAD.jl.

The only argument is a String containing the path to
nomad.3.9.1 folder.

"""
function init(path_to_nomad::String)
	@info "loading NOMAD libraries"
	nomad_libs_call(path_to_nomad)
	create_Evaluator_class()
	create_Extended_Poll_class()
	create_Cresult_class()
	create_cxx_runner()
end

"""

	nomad_libs_call(".../nomad.3.9.1")

load sgtelib and nomad libraries needed to run NOMAD.
Also include all headers to access them via Cxx commands.

"""
function nomad_libs_call(path_to_nomad)

	try
		Libdl.dlopen(path_to_nomad * "/lib/libnomad.so", Libdl.RTLD_GLOBAL)
	catch e
		@warn "NOMAD.jl error : initialization failed, cannot access NOMAD libraries, first need to build them"
		throw(e)
	end


	try
		addHeaderDir(joinpath(path_to_nomad,"src"))
		addHeaderDir(joinpath(path_to_nomad,"ext/sgtelib/src"))
		cxxinclude("nomad.hpp")
	catch e
		@warn "NOMAD.jl error : initialization failed, headers folder cannot be found in NOMAD files"
		throw(e)
	end
end

"""

	create_Evaluator_class()

Create a C++-class "Wrap_Evaluator" that inherits from
NOMAD::Evaluator.

"""
function create_Evaluator_class()

	#=
	The method eval_x is called by NOMAD to evaluate the
	values of objective functions and constraints for a
	given state. The first attribute evalwrap of the class
	is a pointer to the julia function that wraps the evaluator
	provided by the user and makes it interpretable by C++.
	This wrapper is called by the method eval_x. This way,
	each instance of the class Wrap_Evaluator is related
	to a given julia evaluator.

	the attribute m is the number of outputs
	(objective functions and constraints).
	=#

    cxx"""
		#include <string>
		#include <limits>
		#include <vector>

		class Wrap_Evaluator : public NOMAD::Evaluator {
		public:

			double * (*evalwrap)(double * input);
			bool sgte;
			int m;

		  Wrap_Evaluator  ( const NOMAD::Parameters & p, double * (*f)(double * input), int output_dim, bool has_sgte) :

		    NOMAD::Evaluator ( p ) {evalwrap=f; m=output_dim; sgte=has_sgte;}

		  ~Wrap_Evaluator ( void ) {evalwrap=nullptr;}

		  bool eval_x ( NOMAD::Eval_Point   & x  ,
				const NOMAD::Double & h_max      ,
				bool                & count_eval   ) const
			{

			int n = x.get_n();

			double c_x[n+2];

			c_x[0] = n; //first coordinate indicates dimension 

			for (int i = 0; i < n; ++i) {
				c_x[i+1]=x[i].value();
			} //next coordinates contain x

			c_x[n+2] = (x.get_eval_type()==NOMAD::SGTE)?1.0:0.0; //last coordinate equal to 1 if call to the surrogate

			double * c_bb_outputs = evalwrap(c_x);
			//table containing result of black box evaluation. two last
			//coordinates are booleans success and count_eval.

			for (int i = 0; i < m; ++i) {
				NOMAD::Double nomad_bb_output = c_bb_outputs[i];
		    	x.set_bb_output  ( i , nomad_bb_output  );
			} //converting C-double returned by evalwrap in NOMAD::Double and
			//inserting them in x as black box outputs

			bool success = false;
			if (c_bb_outputs[m]==1.0) {
				success=true;
			}

			count_eval = false;
			if (c_bb_outputs[m+1]==1.0) {
				count_eval=true;
			}
			//success and count_eval returned by evalwrap are actually doubles and need
			//to be converted to booleans

			delete[] c_bb_outputs;

		    return success;
		}

		};"""
end

"""

	create_Extended_Poll_class()

Create a C++-class "Wrap_Extended_Poll" that inherits from
NOMAD::Extendeed_Poll. This class will only be used for problems
with categorical variables.

"""
function create_Extended_Poll_class()

	#=
	The method construct_extended_points is called by NOMAD to
	provide extended poll points for evaluation. The first
	attribute extpollwrap of the class is a pointer to the julia
	function that wraps the extended poll provided by the user
	and makes it interpretable by C++.
	This wrapper is called by the method construct_extended_points.
	This way,each instance of the class Wrap_Extended_Poll is
	related to a given julia extended poll.

	the attribute signatures is a vector containing the signatures
	provided by the user.
	=#

	cxx"""class Wrap_Extended_Poll : public NOMAD::Extended_Poll
		{

		private:

			double * (*evalwrap)(double * input);

			// signatures
			std::vector<NOMAD::Signature *> s;


		public:

			// constructor:
			Wrap_Extended_Poll ( NOMAD::Parameters & p,  double * (*f)(double * input) , std::vector<NOMAD::Signature *> sign) :

			NOMAD::Extended_Poll ( p ) {s=sign;evalwrap=f;}

			// destructor:
			~Wrap_Extended_Poll ( void ) {

				for (int i = 0; i < s.size(); ++i) {
					delete s[i];
				}

			 }

			// construct the extended poll points:
			void construct_extended_points ( const NOMAD::Eval_Point & x ) {

				int n = x.get_n();

				double c_x[n+1];

				c_x[0] = -n; //negative first coordinate indicates that extended poll points are needed from evalwrap

				for (int i = 0; i < n; ++i) {
					c_x[i+1]=x[i].value();
				}

				double * c_poll_points = evalwrap(c_x);
				//first coordinate of c_poll_points is the number of extended poll points.
				//Then, extended poll points are all concatenated in this
				//double[], each one preceded by the index of its signature.

				int num_pp = c_poll_points[0]; //number of extended poll points

				int index = 1;

				for (int i = 0; i < num_pp; ++i) {
					int sign_index = static_cast<int> ( c_poll_points[index] ); //read signature index
					NOMAD::Signature * pp_sign = s[sign_index]; //load signature in question
					int npp = pp_sign->get_n(); //dimension of poll point
					NOMAD::Point pp (npp);
					for (int j = 0; j < npp; ++j) {
						pp[j] = c_poll_points[index+1+j];
					}
					add_extended_poll_point ( pp , *pp_sign );
					index += npp+1;
				} //Extracting extended poll points from double[] returned by extpollwrap

			}

		};"""

end

"""
	create_cxx_runner()

Create a C++ function cpp_runner that launches NOMAD
optimization process.

"""
function create_cxx_runner()

	#=
	This C++ function takes as arguments a C-instance NOMAD::Parameters
	along with a void pointer to the julia function wrapping the
	evaluator provided by the user (there is another pointer to the
	extended poll in case categorical variables are used). cpp_runner first
	converts the pointer to the appropriate type and then
	constructs an instance of the class Wrap_Evaluator (proceeds the
	same way for Wrap_Extended_Poll if categorical variables are used).
	Afterwards, a MADS instance is run, taking as arguments the Wrap_Evaluator
	and the NOMAD::Parameters instance. At the end, results are extracted
	from the MADS instance.
	=#

    cxx"""
		#include <iostream>
		#include <string>
		#include <list>

		Cresult cpp_runner(NOMAD::Parameters * p,
					NOMAD::Display out,
					int m, //dimension of the output
					void* f_ptr, //pointer to julia evaluator wrapper
					bool has_stat_avg_,
					bool has_stat_sum_,
					bool has_sgte_,
					bool has_extpoll_,
					std::vector<NOMAD::Signature *> signatures
					) {

		  Cresult res;

		  try {

			typedef double * (*fptr)(double * input);

			//conversion from void pointer to appropriate type (ex_ptr is empty without categorical variables)
			fptr f_fun_ptr = reinterpret_cast<fptr>(f_ptr);

		    // custom evaluator creation
		    Wrap_Evaluator ev   ( *p , f_fun_ptr, m, has_sgte_);

			// custom extended poll creation (used only if there are categorical variables)
			Wrap_Extended_Poll ep ( *p , f_fun_ptr, signatures);

			NOMAD::Mads mads ( *p , &ev , &ep , NULL, NULL );

		    // algorithm creation and execution

			mads.run();

			//saving results
			const NOMAD::Eval_Point* bf_ptr = mads.get_best_feasible();
			const NOMAD::Eval_Point* bi_ptr = mads.get_best_infeasible();
			res.set_eval_points(bf_ptr,bi_ptr,m);
			NOMAD::Stats stats;
			stats = mads.get_stats();
			res.bb_eval = stats.get_bb_eval();
			if (has_stat_avg_) {res.stat_avg = (stats.get_stat_avg()).value();}
			if (has_stat_sum_) {res.stat_sum = (stats.get_stat_sum()).value();}
			res.seed = p->get_seed();

			mads.reset();

			res.success = true;

			delete p;

		  }
		  catch ( exception & e ) {
		    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
		  }

		  NOMAD::Slave::stop_slaves ( out );
		  NOMAD::end();

		  return res;
		}
    """
end

"""

	create_Cresult_class()

Create C++ class that store results from simulation.

"""
function create_Cresult_class()
    cxx"""
		class Cresult {
		public:

			//No const NOMAD::Eval_point pointer in Cresult because GC sometimes erase their content

			std::vector<double> bf;
			std::vector<double> bbo_bf;
			std::vector<double> bi;
			std::vector<double> bbo_bi;
			int bb_eval;
			double stat_avg;
			double stat_sum;
			bool success;
			bool has_feasible;
			bool has_infeasible;
			int seed;

			Cresult(){success=false;}

			void set_eval_points(const NOMAD::Eval_Point* bf_ptr,const NOMAD::Eval_Point* bi_ptr,int m){

				has_feasible = (bf_ptr != NULL);

				if (has_feasible) {
					int n_bf = bf_ptr->get_n();
					for (int i = 0; i < n_bf; ++i) {
						bf.push_back(bf_ptr->value(i));
					}
					for (int i = 0; i < m; ++i) {
						bbo_bf.push_back((bf_ptr->get_bb_outputs())[i].value());
					}
				}

				has_infeasible = (bi_ptr != NULL);

				if (has_infeasible) {
					int n_bi = bi_ptr->get_n();
					for (int i = 0; i < n_bi; ++i) {
						bi.push_back(bi_ptr->value(i));
					}
					for (int i = 0; i < m; ++i) {
						bbo_bi.push_back((bi_ptr->get_bb_outputs())[i].value());
					}
				}

			}


		};
	"""
end
