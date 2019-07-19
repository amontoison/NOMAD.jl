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

Create a Cxx-class "Wrap_Evaluator" that inherits from
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

	the attribute n is the dimension of the problem and m
	is the number of outputs (objective functions and
	constraints).
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

			NOMAD::Signature * s = x.get_signature();
			const std::vector<NOMAD::bb_input_type> & it = s->get_input_types();

			int n = x.get_n();

			double c_x[n+2];

			c_x[0] = n;

			for (int i = 1; i <= n; ++i) {
				c_x[i]=x[i-1].value();
			} //first converting our NOMAD::Eval_Point to a double[]

			c_x[n+1] = (x.get_eval_type()==NOMAD::SGTE)?1.0:0.0;

			double * c_bb_outputs = evalwrap(c_x);

			for (int i = 0; i < m; ++i) {
				NOMAD::Double nomad_bb_output = c_bb_outputs[i];
		    	x.set_bb_output  ( i , nomad_bb_output  );
			} //converting C-double returned by evalwrap in NOMAD::Double that
			//are inserted in x as black box outputs

			bool success = false;
			if (c_bb_outputs[m]==1.0) {
				success=true;
			}

			count_eval = false;
			if (c_bb_outputs[m+1]==1.0) {
				count_eval=true;
			}
			//count_eval returned by evalwrap is actually a double and needs
			//to be converted to a boolean

			delete[] c_bb_outputs;

		    return success;
			//the call to eval_x has succeded
		}

		};"""
end

function create_Extended_Poll_class()

	cxx"""class Wrap_Extended_Poll : public NOMAD::Extended_Poll
		{

		private:

			// signatures
			std::vector<NOMAD::Signature *> s;

			double * (*extpollwrap)(double * input);

		public:

			// constructor:
			Wrap_Extended_Poll ( NOMAD::Parameters & p,  double * (*g)(double * input) , std::vector<NOMAD::Signature *> sign) :

			NOMAD::Extended_Poll ( p ) {s=sign;extpollwrap=g;}

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

				c_x[0] = n;

				for (int i = 1; i <= n; ++i) {
					c_x[i]=x[i-1].value();
				} //first converting our NOMAD::Eval_Point to a double[]

				double * c_poll_points = extpollwrap(c_x);
				//first coordinate is the number of extended poll points
				//then extended poll points are all concatenated in this
				//double[], each one preceded by the index of its signature.

				int num_pp = c_poll_points[0]; //number of extended poll points

				int index = 1;

				for (int i = 0; i < num_pp; ++i) {
					int sign_index = static_cast<int> ( c_poll_points[index] );
					NOMAD::Signature * pp_sign = s[sign_index];
					int npp = pp_sign->get_n(); //dimension of poll point
					NOMAD::Point pp (npp);
					for (int j = 0; j < npp; ++j) {
						pp[j] = c_poll_points[index+1+j];
					}
					add_extended_poll_point ( pp , *pp_sign );
					index += npp+1;
				} //Extracting extended poll points from double[] returned by extendwrap

			}

		};"""

end

"""
	create_cxx_runner()

Create a C++ function cpp_main that launches NOMAD
optimization process.

"""
function create_cxx_runner()

	#=
	This C++ function takes as arguments the settings of the
	optimization (dimension, output types, display options,
	bounds, etc.) along with a void pointer to the julia
	function that wraps the evaluator provided by the user.
	cpp_main first create an instance of the C++ class
	Paramaters and feed it with the optimization settings.
	Then a Wrap_Evaluator is constructed from this Parameters
	instance and from the pointer to the evaluator wrapper.
	Mads is then run, taking as arguments the Wrap_Evaluator
	and Parameters instances.
	=#

    cxx"""
		#include <iostream>
		#include <string>
		#include <list>

		Cresult cpp_runner(NOMAD::Parameters * p,
					NOMAD::Display out,
					int m,
					void* f_ptr,
					void* ex_ptr,
					bool has_stat_avg_,
					bool has_stat_sum_,
					bool has_sgte_,
					bool has_extpoll_,
					std::vector<NOMAD::Signature *> signatures
					) { //le C-main prend en entrée les attributs de l'instance julia parameters


			//Attention l'utilisation des std::string peut entrainer une erreur selon la version du compilateur qui a été utilisé pour générer les librairies NOMAD

		  Cresult res;

		  try {

			//conversion from void pointer to appropriate pointer
			typedef double * (*fptr)(double * input);
			fptr f_fun_ptr = reinterpret_cast<fptr>(f_ptr);

		    // custom evaluator creation
		    Wrap_Evaluator ev   ( *p , f_fun_ptr, m, has_sgte_);

			Wrap_Extended_Poll * wrap_ep_ptr=NULL;

			fptr ex_fun_ptr = reinterpret_cast<fptr>(ex_ptr);
			Wrap_Extended_Poll ep ( *p , ex_fun_ptr, signatures);
			wrap_ep_ptr = &ep;

			NOMAD::Mads mads ( *p , &ev , wrap_ep_ptr , NULL, NULL );

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
