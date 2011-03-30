/*
*  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
*  Copyright (C) 2000 Bruno Levy
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy
*
*     levy@loria.fr
*
*     ISA Project
*     LORIA, INRIA Lorraine, 
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX 
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs. 
*/

#ifndef __OGF_MATH_SYMBOLIC_STENCIL__
#define __OGF_MATH_SYMBOLIC_STENCIL__

#include <OGF/math/common/common.h>
#include <OGF/math/symbolic/symbolic.h>
#include <vector>

namespace OGF {

	/**
	* Abstract class for stencils.
	*/
	class MATH_API Stencil {
	public:
		Stencil(
			int nb_variables, int nb_parameters
			) : nb_variables_(nb_variables), nb_parameters_(nb_parameters) {
		}
		virtual ~Stencil() ;
		int nb_variables() const  { return nb_variables_ ;  }
		int nb_parameters() const { return nb_parameters_ ; }
		virtual double f(const Symbolic::Context& args) = 0 ;
		virtual double g(int i, const Symbolic::Context& args) = 0 ;
		virtual double G(int i, int j, const Symbolic::Context& args) = 0 ;
	private:
		int nb_variables_ ;
		int nb_parameters_ ;
	} ;

	namespace Symbolic {

		typedef std::vector<Expression> Gradient ;
		typedef std::vector<std::vector<Expression> > Hessian ;

		/**
		* A Stencil doing formal derivation to compute
		* the gradient and Hessian. 
		*/
		class MATH_API Stencil : public ::OGF::Stencil {
		public:
			Stencil(const Expression& f, bool use_hessian=false) ;

			virtual double f(const Context& args) ;
			virtual double g(int i, const Context& args) ;
			virtual double G(int i, int j, const Context& args) ;

			Expression& f() { return f_ ; }
			Gradient& g() { return gradient_ ; }
			Expression& g(int i) { return gradient_[i] ; }
			Hessian& G() { return hessian_; }
			Expression& G(int i, int j) { return hessian_[ogf_max(i,j)][ogf_min(i,j)]; }
			void print(std::ostream& out) ;
		private:
			Expression f_ ;
			Gradient gradient_ ;
			Hessian  hessian_ ;
			bool use_hessian_;
		} ;
	}

	//________________________________________________________________________________________

	class MATH_API StencilInstance {
	public:
		StencilInstance() : stencil_(nil), global_indices_(nil), parameters_(nil) { }
		/*
		StencilInstance(const StencilInstance& rhs) { 
			stencil_=rhs.stencil_;global_indices_=rhs.global_indices_;
			parameters_ = rhs.parameters_;
		}
		StencilInstance& operator=(const StencilInstance& rhs){
			stencil_=rhs.stencil_;global_indices_=rhs.global_indices_;
			parameters_ = rhs.parameters_;return *this;}*/
		~StencilInstance() ;

		int nb_variables() const { return stencil_->nb_variables() ; }
		int nb_parameters() const { return stencil_->nb_parameters() ; }
		int& global_variable_index(int i) { 
			assert(i >= 0 && i < stencil_->nb_variables()) ; 
			return global_indices_[i] ; 
		}
		double& parameter(int i) {
			assert(i >= 0 && i < stencil_->nb_parameters()) ; 
			return parameters_[i] ;
		}
		Stencil* stencil() { return stencil_ ; }
public:
		void bind(Stencil* stencil) ;
		bool is_initialized() ;
		
	private:
		Stencil* stencil_ ;
		int* global_indices_ ;
		double* parameters_ ;
		
	} ;
}

#endif