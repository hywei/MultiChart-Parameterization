#include <OGF/math/symbolic/stencil.h>

namespace OGF {

	//__________________________________________________________________________________________________

	Stencil::~Stencil() {
	}

	//__________________________________________________________________________________________

	namespace Symbolic {

		Stencil::Stencil(
			const Expression& f, bool use_hessian
			) : ::OGF::Stencil(f->max_variable_index() + 1, f->max_parameter_index()+1),
			f_(f), 
			gradient_(nb_variables())
		{
			bool ok = true ;
			use_hessian_ = use_hessian;

			for(int i=0; i<nb_variables(); i++) {
				if(!f_->depends_on_variable(i)) {
					std::cerr << "variable " << i << " not used in stencil" << std::endl ;
					ok = false ;
				}
			}


			for(int i=0; i<nb_parameters(); i++) {
				if(!f_->depends_on_parameter(i)) {
					std::cerr << "parameter " << i << " not used in stencil" << std::endl ;
					ok = false ;
				}
			}

			assert(ok) ;

			for(int i=0; i<nb_variables(); i++) {
				gradient_[i] = der(f_, i) ;
			}

			if (use_hessian) {
				hessian_.resize(nb_variables());

				for(int i=0; i<nb_variables(); i++) {
					hessian_[i].resize(i+1);
					for(int j=0; j<=i; j++) {
						hessian_[i][j] = der(gradient_[i], j) ;
					}
				}
			}
		}

		void Stencil::print(std::ostream& out) {
			out << "f=" << f_ << std::endl ;
			out << std::endl ;
			{
				for(int i=0; i<nb_variables(); i++) {
					out << "g" << i << "=" << gradient_[i] << std::endl ;
				}
			}

			if (use_hessian_){
				out << std::endl ;
				{
					for(int i=0; i<nb_variables(); i++) {
						for(int j=0; j<=i; j++) {
							out << "G" << i << "," << j << "=" << hessian_[i][j] << std::endl ;
						}
					}
				}
			}
			out << std::endl ;
		}

		double Stencil::f(const Symbolic::Context& args) {
			return f()->eval(args) ;
		}

		double Stencil::g(int i, const Symbolic::Context& args) {
			return g(i)->eval(args) ;
		}

		double Stencil::G(int i, int j, const Symbolic::Context& args) {
			return G(i,j)->eval(args) ;
		}
	}

	//__________________________________________________________________________________________________

	void StencilInstance::bind(Stencil* stencil) {
		assert(stencil_ == nil) ;
		assert(global_indices_ == nil) ;
		assert(parameters_ == nil) ;
		stencil_ = stencil ;
		global_indices_ = new int[stencil_->nb_variables()] ;
		parameters_     = new double[stencil_->nb_parameters()] ;
		{for(int i=0; i<stencil_->nb_variables(); i++) {
			global_indices_[i] = -1 ;
		}}
		{for(int i=0; i<stencil_->nb_parameters(); i++) {
			parameters_[i] = 0 ;
		}}
	}

	bool StencilInstance::is_initialized() {
		if(stencil_ == nil) { return false ; }
		for(int i=0; i<stencil_->nb_variables(); i++) {
			if(global_indices_[i] == -1) {
				std::cerr << "global index:" << i << " : unitialized" << std::endl ;
				return false ;
			}
		}
		return true ;
	}

	StencilInstance::~StencilInstance() { 
		if (global_indices_ != nil){
			delete[] global_indices_ ; 
			global_indices_ = nil ; 
		}
		if (parameters_ != nil){
			delete[] parameters_ ;
			parameters_ = nil ;
		}
		stencil_ = nil ;
	}
}