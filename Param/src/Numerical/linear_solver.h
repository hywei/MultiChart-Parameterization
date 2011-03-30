//
// This class is the implementation of a row-dominant sparse matrix.
// Primarily, the class is designed for sparse mesh structure.
//
//////////////////////////////////////////////////////////////////////

#ifndef LINEAR_SOLVER_2_H
#define LINEAR_SOLVER_2_H

#include "MeshSparseMatrix.h"
#include "solver.h"
#ifdef WIN32
#include <hj_3rd/hjlib/sparse_old/sparse.h>
#else
#include <hj_3rd/hjlib/sparse/sparse.h>
#endif

class LinearSolver: public Solver
{
public:
	LinearSolver(int nb_variables);
    virtual ~LinearSolver();

public:
	// __________________ Construction _____________________
	void begin_equation() ;

	void begin_row() ;
	void set_right_hand_side(double b_) ;
	void add_coefficient(int index_, double a_) ;
	void end_row() ;
	void end_equation() ;

	//
	void solve();
	void solve_from_file();

	//
	void factorize();
	void renew_right_b(std::vector<double>& right_b_vec);
	void set_equation_div_flag();

	void equations_value(std::vector<double>& var_val_vec);
	size_t get_equation_size(){return m_equation_vec.size();}

	void is_printf_info(bool is_){m_is_printf_info=is_;}
	void print_to_file(std::vector<double>& var_val_vec, std::string filename);
	void write_to_file(std::string ata_filename, std::string atb_filename);

private:

	void set_solve_matrix();
	void set_solve_b();
	void update_variables();

	bool is_free(int id)   { return (id < nb_free_variables_) ;  }
	bool is_locked(int id) { return (id >= nb_free_variables_) ; }

	void print_f(std::vector<double>& xc_);
	void print_equation_value(std::vector<double>& function_vec);

private:
	// User representation
	int nb_free_variables_ ;
	int nb_locked_variables_ ;

    std::vector<std::vector<std::pair<int, double> > > m_equation_vec;
    std::vector<std::pair<int, double> > m_current_equ;
    std::vector<double> m_right_b_vec;
    std::vector<double> m_xc_;

private:
	//std::auto_ptr<hj::sparse::solver> m_solver_;
	CMeshSparseMatrix m_solve_matrix_AT_;
    std::vector<double> m_solve_b_vec;

private:
	bool factorize_state;

private:
    std::vector<size_t> m_equ_div_flag_vec;
	bool m_is_printf_info;

    std::vector<std::pair<int, double> > m_tmp_current_equ;
    std::vector<std::vector<std::pair<int, double> > > m_tmp_equation_vec;
};

#endif
