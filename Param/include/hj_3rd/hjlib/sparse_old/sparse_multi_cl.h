#ifndef _SPARSEMULTICL_H_
#define _SPARSEMULTICL_H_

#include "sparse.h"

// AT must be pre-allocated
void HJ_SPARSE_API trans(const hj::sparse::spm_csc<double> &A, hj::sparse::spm_csc<double> &AT);

inline void transx(const hj::sparse::spm_csc<double> &A, hj::sparse::spm_csc<double> &AT)
{
	AT.resize(A.size(2), A.size(1), hj::sparse::nnz(A));
	trans(A, AT);
}

struct mm_rtn
{
	int rows, cols, nnz;
	void * ctx;
};

mm_rtn HJ_SPARSE_API spm_dmm(bool is_A_trans, const hj::sparse::spm_csc<double> &A,
			   bool is_B_trans, const hj::sparse::spm_csc<double> &B, mm_rtn *rtn = 0);

bool HJ_SPARSE_API convert(const mm_rtn &rtn, hj::sparse::spm_csc<double> &M);

void HJ_SPARSE_API mm_rtn_destroy(mm_rtn &rtn);

inline bool spm_dmm(bool is_A_trans, const hj::sparse::spm_csc<double> &A,
			   bool is_B_trans, const hj::sparse::spm_csc<double> &B,
			   hj::sparse::spm_csc<double> &C, mm_rtn *rtn = 0)
{
	mm_rtn rtn1 = spm_dmm(is_A_trans, A, is_B_trans, B, rtn);
	if(!rtn1.ctx)
		return false;
	C.resize(rtn1.rows, rtn1.cols, rtn1.nnz);
	bool rtn2 = convert(rtn1, C);
	if(rtn == 0)
		mm_rtn_destroy(rtn1);
	return rtn2;
}

#endif
