#ifndef HJ_SPARSE_DEPRECATED_H_
#define HJ_SPARSE_DEPRECATED_H_

#include "config.h"

#include "format.h"
#include "solver.h"
#include "operation.h"

#include <memory>

#ifdef __GNUG__
#define HJ_SPARSE_DEPRECATED __attribute__((deprecated))
#endif

#ifdef _MSC_VER
#define HJ_SPARSE_DEPRECATED __declspec(deprecated)
#endif

namespace hj { namespace sparse {

//! @brief csc
template <typename VAL_TYPE= double, typename INT_TYPE = idx_type>
class HJ_SPARSE_DEPRECATED spm_csc : public csc<VAL_TYPE, INT_TYPE>
{
public:
 spm_csc()
:csc<VAL_TYPE, INT_TYPE>(0, 0, 0) {
  }
	spm_csc(INT_TYPE rows, INT_TYPE cols, INT_TYPE nnz = 0)
		:csc<VAL_TYPE, INT_TYPE>(rows, cols, nnz) {
	}
};

class solver
{
public:
	/**
	   @param id "umfpack" or "cholmod", default 0 is current best
	   implementation
	*/
	HJ_SPARSE_DEPRECATED static solver *create(const csc<double, int> &A, const char *id = 0) {
		std::auto_ptr<direct_solver_A> slv(direct_solver_A::create(A, id));
		if(slv.get())
			return new solver(slv.release());
		return 0;
	}
	virtual bool solve(const double *b, double *x, int nrhs = 1) {
		return !slv_->solve(b, x, nrhs);
	}
private:
	solver(direct_solver_A *slv)
		:slv_(slv) {
	}
	std::auto_ptr<direct_solver_A> slv_;
};

template <typename T1, typename INT_TYPE1,
		  typename T2, typename INT_TYPE2,
		  typename T3, typename INT_TYPE3>
HJ_SPARSE_DEPRECATED
inline csc<T3, INT_TYPE3> &spm_mpm(const csc<T1, INT_TYPE1> &A, const csc<T2, INT_TYPE2> &B, csc<T3, INT_TYPE3> &C)
{
	csc_by_vm<T3, INT_TYPE3, std_map> vmC(A.size(1), A.size(2)); //! NOTICE: nz per col
	acc(A, vmC); // equal to convert A to vmC
	acc(B, vmC);
	convert(vmC, C);
	return C;
}

template <typename T, typename INT_TYPE, typename M1, typename M2>
HJ_SPARSE_DEPRECATED
M2 &mm(bool transA, const csc<T, INT_TYPE> &A,
	   bool transB, const zjucad::matrix::matrix_expression<M1> &B,
	   bool transC, zjucad::matrix::matrix_expression<M2> &C)
{
	using namespace zjucad::matrix;
#define spm_mm_MULT(op) \
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) { for(; vi < A.ptr_[ci+1]; ++vi) { \
		op;\
	}}

	INT_TYPE vi = A.ptr_[0];
	if(transA) {
		if(transB) {
			if(transC) {
				spm_mm_MULT(C()(colon(), ci) += A.val_[vi]*B()(colon(), A.idx_[vi]));
			}
			else {
				spm_mm_MULT(C()(ci, colon()) += A.val_[vi]*trans(B()(colon(), A.idx_[vi])));
			}
		}
		else {
			if(transC) {
				spm_mm_MULT(C()(colon(), ci) += A.val_[vi]*trans(B()(A.idx_[vi], colon())));
			}
			else {
				spm_mm_MULT(C()(ci, colon()) += A.val_[vi]*B()(A.idx_[vi], colon()));
			}
		}
	}
	else {
		if(transB) {
			if(transC) {
				spm_mm_MULT(C()(colon(), A.idx_[vi]) += A.val_[vi]*B()(colon(), ci));
			}
			else {
				spm_mm_MULT(C()(A.idx_[vi], colon()) += A.val_[vi]*trans(B()(colon(), ci)));
			}
		}
		else {
			if(transC) {
				spm_mm_MULT(C()(colon(), A.idx_[vi]) += A.val_[vi]*trans(B()(ci, colon())));
			}
			else {
				spm_mm_MULT(C()(A.idx_[vi], colon()) += A.val_[vi]*B()(ci, colon()));
			}
		}
	}
#undef spm_mm_MULT
	return C();
}

template <typename T, typename INT_TYPE, typename M1, typename M2>
HJ_SPARSE_DEPRECATED
M2 &spm_mm(bool transA, const csc<T, INT_TYPE> &A,
		   bool transB, const zjucad::matrix::matrix_expression<M1> &B,
		   bool transC, zjucad::matrix::matrix_expression<M2> &C)
{
	return mm(transA, A, transB, B, transC, C);
}

template <typename T, typename INT_TYPE, typename V1, typename V2>
HJ_SPARSE_DEPRECATED
V2 &spm_mxv(const csc<T, INT_TYPE> &A, const V1 &b, V2 &c)
{
	INT_TYPE vi = A.ptr_[0];
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
		for(; vi < A.ptr_[ci+1]; ++vi) {
			c[A.idx_[vi]] += A.val_[vi]*b[ci];
		}
	}
	return c;
}

template <typename T, typename INT_TYPE, typename V1, typename V2>
HJ_SPARSE_DEPRECATED
V2 &spm_mTxv(const csc<T, INT_TYPE> &A, const V1 &b, V2 &c)
{
	INT_TYPE vi = A.ptr_[0];
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
		for(; vi < A.ptr_[ci+1]; ++vi) {
			c[ci] += A.val_[vi]*b[A.idx_[vi]];
		}
	}
	return c;
}

}}

/**
  rtn can be used for caching the patten
  */
HJ_SPARSE_DEPRECATED
inline void spm_dmm(bool is_A_trans, const hj::sparse::csc<double> &A,
			   bool is_B_trans, const hj::sparse::csc<double> &B,
					hj::sparse::csc<double> &C)
{
	hj::sparse::MM<>(is_A_trans, A, is_B_trans, B, C);
}

#undef HJ_SPARSE_DEPRECATED

#endif
