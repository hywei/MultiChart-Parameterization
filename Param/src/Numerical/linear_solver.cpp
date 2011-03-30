#include "linear_solver.h"
#include "../Common/stopwatch.h"
#include "../Numerical/MatrixConverter.h"

#ifdef WIN32
#include <hj_3rd/hjlib/sparse_old/sparse_multi_cl.h>
#else
#include <hj_3rd/hjlib/sparse/sparse_multi_cl.h>
#endif

#include <assert.h>
#include <float.h>
#include <set>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <memory>
using namespace std;
//#include <process.h>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
LinearSolver::LinearSolver(int nb_variables) 
: Solver(nb_variables) 
{
	factorize_state = false;
	m_is_printf_info = true;
}
LinearSolver::~LinearSolver()
{
}
//////////////////////////////////////////////////////////////////////
// public methods
//////////////////////////////////////////////////////////////////////
void LinearSolver::begin_equation() 
{
	// check if referred variable is also setted value
	int valued_var_num;
	do 
	{
		valued_var_num = 0;
		for(int i=0; i<nb_variables_; i++) {
			
			if (variable_[i].is_referred() && variable_[i].is_valued()) {
				int ref_var_index = variable_[i].ref_var_index();
				if (variable_[ref_var_index].is_valued()) {
					double val_ = variable_[ref_var_index].value();
					if (val_ != variable_[i].value()) {
						printf("%d and %d variable value are conflicted.\n", i, ref_var_index);
					}
				} else {
					variable_[ref_var_index].lock();
					variable_[ref_var_index].set_value(variable_[i].value());
					valued_var_num++;
				}
			}
		}
	} while(valued_var_num > 0);

	for(int i=0; i<nb_variables_; i++) {
		
		if (variable_[i].is_referred() && variable_[i].is_valued()) {
			variable_[i].set_ref_variable_inex(-1);
		}
	}

	//printf("nb_variables_ = %d\n", nb_variables_);

	//
	nb_free_variables_ = 0 ;
	nb_locked_variables_ = 0 ;
	int cur_index = 0 ;
	for(int i=0; i<nb_variables_; i++) {
		
		if(!variable_[i].is_locked()) {
			set_variable_index(variable_[i], cur_index) ;
			nb_free_variables_++ ;
			cur_index++ ;
		}
	}
	//printf("nb_free_variables_ = %d, \t nb_locked_variables_ = %d\n", nb_free_variables_, nb_locked_variables_);
	
	//
	for(int i=0; i<nb_variables_; i++) {
	
		if(variable_[i].is_valued()) {
			
			set_variable_index(variable_[i], cur_index) ;
			cur_index++ ;
			nb_locked_variables_++ ;
		}
	}
	//printf("nb_free_variables_ = %d, \t nb_locked_variables_ = %d\n", nb_free_variables_, nb_locked_variables_);
	for(int i=0; i<nb_variables_; i++) {
		
		if (variable_[i].is_referred()) {
			int ref_var_index = variable_[i].ref_var_index();
			if (variable_[ref_var_index].is_valued()) {
				
				variable_[i].set_value(variable_[ref_var_index].value());
				set_variable_index(variable_[i], cur_index) ;
				cur_index++ ;
				nb_locked_variables_++ ;
			} else {
				set_variable_index(variable_[i], variable_[ref_var_index].index());
				variable_[i].unlock();
			}
		}
	}

	m_equ_div_flag_vec.push_back(0);
}
void LinearSolver::begin_row() 
{
	m_current_equ.clear();
	//m_tmp_current_equ.clear();
}
void LinearSolver::set_right_hand_side(double b_) 
{
	m_right_b_vec.push_back(b_);
}
void LinearSolver::add_coefficient(int index_, double a_) 
{
	int internal_id = variable_[index_].index() ;
	m_current_equ.push_back(make_pair(internal_id, a_));

	//m_tmp_current_equ.push_back(make_pair(index_, a_));
}
void LinearSolver::end_row() 
{
	m_equation_vec.push_back(m_current_equ);

// 	for(size_t k=0; k<m_current_equ.size(); ++k)
// 	{
// 		std::cout << m_current_equ[k].second <<" ";
// 	}
// 	std::cout << std::endl;
	//m_tmp_equation_vec.push_back(m_tmp_current_equ);
}

void LinearSolver::end_equation()
{
	m_xc_.resize(nb_variables_) ;
	for(int i=0; i<nb_variables_; i++) {
	if(variable_[i].index() <0 || variable_[i].index() > nb_variables_*2) 
			printf("i = %d, idx = %d, value = %d\n", i, variable_[i].index(), variable_[i].value());

		m_xc_[variable_[i].index()] = variable_[i].value() ;
	}
	m_equ_div_flag_vec.push_back(m_equation_vec.size());

	//
/*	printf("set solve b!\n");*/
	set_solve_b();
/*	printf("set solve matrix!\n");*/
	set_solve_matrix();
}
void  LinearSolver::solve_from_file()
{
	// write_to_file("c:\\solver\\ata.tmp", "c:\\solver\\atb.tmp");

	// printf("ouput to file.\n");

	// /*
	// PROCESS_INFORMATION piProcInfo; 
	// STARTUPINFO siStartInfo;

	// // Set up members of STARTUPINFO structure.
	// siStartInfo.cb = sizeof(STARTUPINFO); 
	// siStartInfo.lpReserved = NULL;
	// siStartInfo.lpReserved2 = NULL; 
	// siStartInfo.cbReserved2 = 0;
	// siStartInfo.lpDesktop = NULL; 
	// siStartInfo.dwFlags = 0;

	// CreateProcess(NULL, (LPWSTR)"c:\\solver\\solver.exe ata.tmp atb.tmp mx.tmp",
	// 	NULL,NULL,0,0,NULL,NULL,&siStartInfo, &piProcInfo);

	// // Wait for the processs to finish
	// DWORD rc = WaitForSingleObject(
	// 	piProcInfo.hProcess, // process handle
	// 	INFINITE); */

	// char *argv[] = {
	// 	"c:\\solver\\solver.exe",
	// 	"c:\\solver\\ata.tmp",
	// 	"c:\\solver\\atb.tmp",
	// 	"c:\\solver\\mx.tmp",
	// 	0
	// };
	// int rst =spawnvp(_P_WAIT, "c:\\solver\\solver.exe", argv);

	// printf("solve = %d\n", rst);

	// printf("read from file.\n");
	// //
	// ifstream ifs("c:\\solver\\mx.tmp", std::ios::binary);
	// if (ifs.fail())
	// {
	// 	printf("can't load mx.\n");
	// }
	// zjucad::matrix::matrix<double> mx_matrix;
	// printf("%d\n", __LINE__);
	// read(ifs, mx_matrix);
	// printf("%d\n", __LINE__);

	// vector<double> m_x_(mx_matrix.size());
	// printf("%d\n", __LINE__);
	// for (int i = 0; i < mx_matrix.size(); i++)
	// {
	// 	m_x_[i] = mx_matrix(i);
	// }
	// printf("%d\n", __LINE__);

	// printf("m_xc size = %d\n", m_xc_.size());
	// printf("m_x_ size = %d\n", m_x_.size());

	// //
	// for(int i=0; i<nb_free_variables_; i++) {
	// 	m_xc_[i] = m_x_[i] ;
	// }

	// //	
	// update_variables();

	// if (m_is_printf_info){
	// 	printf("after solving: \n");
	// 	//print_f(m_x_);
	//}
}
void LinearSolver::solve()
{

	//factorize();
	vector<double> at_b_vec;
	vector<double> m_x_(m_solve_matrix_AT_.GetRowNum());
	m_solve_matrix_AT_.MultiplyVector(m_solve_b_vec, at_b_vec);

	//
	std::auto_ptr<hj::sparse::solver> m_solver_;

	// H = JT * J
	SystemStopwatch watch_;
	{
        static hj::sparse::csc<double,int> spm_ATA;

		hj::sparse::csc<double> spm_A;
		CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_A, m_solve_matrix_AT_);

		hj::sparse::MM<>(false, spm_A, true, spm_A, spm_ATA);

		m_solve_matrix_AT_.ClearData();
		m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));
	}
	//
	if (m_is_printf_info){
		watch_.print_elapsed_time();
	}

	//
	m_equation_vec.clear();
	if(!m_solver_.get()) {
		printf("factorize failed.\n");
	}
	m_solve_b_vec.clear();
	m_right_b_vec.clear();
	bool su = m_solver_->solve(&at_b_vec[0], &m_x_[0]);
	
	if (!su)
	{
		printf("solve x failed.!!!\n");
	}

	//
	for(int i=0; i<nb_free_variables_; i++) {
		m_xc_[i] = m_x_[i] ;
	}

	//	
	update_variables();

	if (m_is_printf_info){
		printf("after solving: \n");
		//print_f(m_x_);
	}
}

void LinearSolver::factorize()
{
	/*
	if (!factorize_state)
	{
		// H = JT * J
		SystemStopwatch watch_;
		hj::sparse::spm_csc<double> spm_ATA;
		{
			hj::sparse::spm_csc<double> spm_A;
			CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_A, m_solve_matrix_AT_);

			spm_dmm(false, spm_A, true, spm_A, spm_ATA);
		}
		m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));

		if (m_is_printf_info){
			watch_.print_elapsed_time();
		}
		
		factorize_state = true;
	}*/
}
void LinearSolver::renew_right_b(vector<double>& right_b_vec)
{
	//m_right_b_vec.assign(right_b_vec.begin(), right_b_vec.end());
 	m_right_b_vec = right_b_vec;
	set_solve_b();
}
void LinearSolver::equations_value(vector<double>& var_val_vec)
{
	vector<double> input_x_(nb_free_variables_);

	//
	for(int i=0; i<nb_variables_; i++) {
		if (is_free(variable_[i].index())) {
			input_x_[variable_[i].index()] = var_val_vec[i] ;
		}
	}

	print_f(input_x_);
}
void LinearSolver::set_equation_div_flag()
{
	m_equ_div_flag_vec.push_back(m_equation_vec.size());
}
//////////////////////////////////////////////////////////////////////
// private methods
//////////////////////////////////////////////////////////////////////
void LinearSolver::set_solve_matrix()
{
	//
	vector<bool> row_valid_flag(m_equation_vec.size());
	fill(row_valid_flag.begin(), row_valid_flag.end(), false);
	int row_size = 0;
	for (size_t i = 0; i < m_equation_vec.size(); i++)
	{
		vector<pair<int, double> >& cur_equ = m_equation_vec[i];

		//
		size_t lock_num = 0;
		for (size_t j = 0 ; j < cur_equ.size(); j++)
		{
			pair<int, double>& var_ = cur_equ[j];
			if (is_locked(var_.first))
			{
				lock_num++;
			}
		}
		if (lock_num < cur_equ.size()) {
			row_size++;
			row_valid_flag[i] = true;
		}
	}
	int col_size = nb_free_variables_;

	//
	CMeshSparseMatrix solve_matrix;
	solve_matrix.SetRowCol(row_size, col_size);

	int row_ = 0;
	for (size_t i = 0; i < m_equation_vec.size(); i++)
	{
		if (row_valid_flag[i])
		{
			vector<pair<int, double> >& cur_equ = m_equation_vec[i];

			//
			for (size_t j = 0 ; j < cur_equ.size(); j++)
			{
				pair<int, double>& var_ = cur_equ[j];
				if (is_free(var_.first)) 
				{
                    // Warning : the following comment by hywei, there may be error
					// if (_isnan(var_.second))
					// {
					// 	printf("is nan.\n");
					// 	size_t aaaaaa = 6;
					//}
					solve_matrix.AddElement(row_, var_.first, var_.second);
				}
			}
			row_++;
		}
	}

	//
	size_t invalid_num = 0;
	vector<size_t> tmp_equ_div_flag_vec(m_equ_div_flag_vec);
	for (size_t i = 1; i < tmp_equ_div_flag_vec.size(); i++)
	{
		size_t start_eqn = tmp_equ_div_flag_vec[i-1];
		size_t end_eqn = tmp_equ_div_flag_vec[i];
		for (size_t j = start_eqn; j < end_eqn; j++)
		{
			if (!row_valid_flag[j])
			{
				invalid_num++;
			}
		}
		m_equ_div_flag_vec[i] = end_eqn - invalid_num;
	}

	//
	solve_matrix.Transpose(m_solve_matrix_AT_);
}
void LinearSolver::set_solve_b()
{
	m_solve_b_vec.clear();
	for (size_t i = 0; i < m_equation_vec.size(); i++)
	{
		vector<pair<int, double> >& cur_equ = m_equation_vec[i];

		if (i == 1251)
		{
			size_t aaaaa = 5;
		}

		//
		double sum_b = m_right_b_vec[i];
		size_t lock_num = 0;
		for (size_t j = 0 ; j < cur_equ.size(); j++)
		{
			pair<int, double>& var_ = cur_equ[j];
			if (is_locked(var_.first))
			{
				sum_b -= m_xc_[var_.first] * var_.second;
				lock_num++;
			}
		}

		if (lock_num < cur_equ.size())
		{
          // Warning: the following codes commented by hywei, there may be error 
			// if (_isnan(sum_b))
			// {
			// 	printf("is nan b.\n");
			// 	size_t aaaaaa = 6;
			// }
			m_solve_b_vec.push_back(sum_b);
			
		}
	}
}
void LinearSolver::update_variables()
{
	for(int i=0; i<nb_variables_; i++) {
		variable_[i].set_value(m_xc_[variable_[i].index()]) ;
	}
}
void LinearSolver::print_f(vector<double>& xc_)
{
	//
	CMeshSparseMatrix solve_matrix_A;
	m_solve_matrix_AT_.Transpose(solve_matrix_A);
	vector<double> function_vector(solve_matrix_A.GetRowNum());
	solve_matrix_A.MultiplyVector(xc_, function_vector);
	for (size_t i = 0; i < function_vector.size(); i++)
	{
		function_vector[i] -= m_solve_b_vec[i];
	}

	print_equation_value(function_vector);
}
void LinearSolver::print_equation_value(vector<double>& function_vec)
{
	double f2_sum_sum = 0;
	for (size_t i = 1; i < m_equ_div_flag_vec.size(); i++)
	{
		size_t start_equ = m_equ_div_flag_vec[i-1];
		size_t end_equ = m_equ_div_flag_vec[i];

		double f2_sum = 0;
		for (size_t j = start_equ; j < end_equ; j++)
		{
			f2_sum += function_vec[j] * function_vec[j];
		}
		f2_sum_sum += f2_sum;
		printf("%d part: %f  ", i-1, f2_sum);
	}
	printf("\nfunction value: %f \n", f2_sum_sum);
}
void LinearSolver::print_to_file(vector<double>& var_val_vec, string filename)
{
	vector<double> input_x_(nb_free_variables_);

	//
	for(int i=0; i<nb_variables_; i++) {
		if (is_free(variable_[i].index())) {
			input_x_[variable_[i].index()] = var_val_vec[i] ;
		}
	}
	CMeshSparseMatrix solve_matrix_A;
	m_solve_matrix_AT_.Transpose(solve_matrix_A);
	vector<double> function_vector(solve_matrix_A.GetRowNum());
	solve_matrix_A.MultiplyVector(input_x_, function_vector);
	for (size_t i = 0; i < function_vector.size(); i++)
	{
		function_vector[i] -= m_solve_b_vec[i];
	}

	size_t end_equ = m_equ_div_flag_vec[1];
	ofstream file(filename.c_str());
	file << end_equ << '\n';

	

	for (size_t i = 0; i < end_equ; i++)
	{
		file << i << ' ' << function_vector[i] << '\n';
	}
	file.close();
}
void LinearSolver::write_to_file(string ata_filename, string atb_filename)
{
	// vector<double> at_b_vec;
	// vector<double> m_x_(m_solve_matrix_AT_.GetRowNum());
	// m_solve_matrix_AT_.MultiplyVector(m_solve_b_vec, at_b_vec);


	// hj::sparse::spm_csc<double> spm_ATA;
	// hj::sparse::spm_csc<double> spm_A;
	// CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_A, m_solve_matrix_AT_);
	// spm_dmm(false, spm_A, true, spm_A, spm_ATA);

	// //
	// zjucad::matrix::matrix<double> atb_(1, (int)at_b_vec.size());
	// for (size_t i = 0; i < at_b_vec.size(); i++)
	// {
	// 	atb_(0, i) = at_b_vec[i];
	// }

	// //
	// ofstream ata_ofs(ata_filename.c_str(), std::ios::binary);
	// //ofstream ata_ofs(ata_filename.c_str());
	// write(ata_ofs, spm_ATA.ptr_);
	// write(ata_ofs, spm_ATA.idx_);
	// write(ata_ofs, spm_ATA.val_);
	// ata_ofs.close();

	// //
	// ofstream atb_ofs(atb_filename.c_str(), std::ios::binary);
	// //ofstream atb_ofs(atb_filename.c_str());
	// write(atb_ofs, atb_);
	// atb_ofs.close();
}
