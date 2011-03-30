#ifndef __OGF_MATH_NUMERIC_SOLVER__
#define __OGF_MATH_NUMERIC_SOLVER__

#include <assert.h>
#include <vector>

class SolverVariable 
{
public:
	SolverVariable() : x_(0), ref_var_index_(-1), index_(-1), 
		               locked_(false), valued_(false) { }

	double value() const { return x_; }
	void set_value(double x_in) { x_ = x_in; valued_=true; /*printf("set index %d valued!\n", index_);*/ }
	void set_ref_variable_inex(int ref_index){ref_var_index_=ref_index;  /*printf("set index %d ref!\n", index_);*/ }

	void lock()   { locked_ = true ; }
	void unlock() { locked_ = false ; }
	bool is_locked() const { return locked_ ; }
	bool is_valued() const {return valued_;}

	bool is_referred() const {return ref_var_index_>=0;}
	int ref_var_index() const {return ref_var_index_;}

	int index() const { return index_ ;}
	void set_index(int index_in) { index_ = index_in ; }

private:
	double x_ ;
	int ref_var_index_;
	int index_ ;
	bool locked_ ;
	bool valued_;

	friend class Solver ;
} ;

class Linear_Constraint
{
public:
	Linear_Constraint() : m_right_b_(0) {}

	//
	void begin_constraint() {m_equation.clear();}
	void set_right_hand_side(double b_) {m_right_b_=b_;}
	void add_coefficient(int index_, double a_) {m_equation.push_back(std::make_pair(index_, a_));}

	//
	double right_b() {return m_right_b_;}
    std::vector<std::pair<int, double> > get_constraint_coef() {return m_equation;}
private:
    std::vector<std::pair<int, double> > m_equation;
	double m_right_b_;

	friend class Solver ;
};

//____________________________________________________________________________

 class Solver {
    public:
        Solver(int nb_variables) ;
        virtual ~Solver() ;

        int nb_variables() const { return nb_variables_ ; }
        
        SolverVariable& variable(int idx) { 
            assert(idx >= 0 && idx < nb_variables_) ;
            return variable_[idx] ;
        }
            
        const SolverVariable& variable(int idx) const {
            assert(idx >= 0 && idx < nb_variables_) ;
            return variable_[idx] ;
        }
        
        void set_variable_index(SolverVariable& var, int index) {
            var.set_index(index) ;
        }

		virtual void solve()=0;

    protected:
        // User representation
        int nb_variables_ ;
        SolverVariable* variable_ ;
        bool ok_ ;
    } ;

//____________________________________________________________________________

#endif
