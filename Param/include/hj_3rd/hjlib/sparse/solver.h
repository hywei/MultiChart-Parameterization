#ifndef HJ_SPARSE_SOLVER_H_
#define HJ_SPARSE_SOLVER_H_

#include "config.h"

extern "C" {

//! sizeof_val = sizeof(float,double), not sizeof(complex<float>, complex<double>)

///! direct solver
HJ_SPARSE_API
void *hj_sparse_direct_solver_A_create(
	unsigned char sizeof_int, 	unsigned char sizeof_val, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *ptr, const void *idx, const void *val, //< link for internal use
	const char *name, //< copy to internal
	void *opts); //< link for internal use

HJ_SPARSE_API
int hj_sparse_direct_solver_A_solve(void *ctx, const void *b, void *x, const size_t nrhs, void *opts);

HJ_SPARSE_API
void hj_sparse_direct_solver_A_destroy(void *ctx);


///! iterative solver
HJ_SPARSE_API
void *hj_sparse_iterative_solver_A_create(
	unsigned char sizeof_int, 	unsigned char sizeof_val, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *ptr, const void *idx, const void *val, //< link for internal use
	const char *name, //< copy to internal
	void *opts); //< link for internal use

HJ_SPARSE_API
int hj_sparse_iterative_solver_A_solve(void *ctx, const void *b, void *x, const size_t nrhs, void *opts);

HJ_SPARSE_API
void hj_sparse_iterative_solver_A_destroy(void *ctx);

struct laspack_opts
{
	char *iter_name_, *precond_name_;
	size_t iter_num_;
	double precision_, relax_;
};

}

#include <complex>

#include "./type_traits.h"

namespace hj { namespace sparse {

class solver_base
{
public:
	virtual ~solver_base(){}
	virtual int solve(const void *b, void *x, size_t nrhs = 1, void *opts = 0) = 0;
protected:
	solver_base(void *ctx):ctx_(ctx){}
	void *ctx_;
};

//! direct solvers

//! @brief solve Ax=b
class direct_solver_A : public solver_base
{
public:
	template <typename VAL_TYPE, typename INT_TYPE>
	static direct_solver_A *create(
		const csc<VAL_TYPE, INT_TYPE> &A, const char *name = 0, void *opts = 0) {
		void *ctx = hj_sparse_direct_solver_A_create(
			sizeof(INT_TYPE), value_type<VAL_TYPE>::size, value_type<VAL_TYPE>::type,
			A.size(1), A.size(2),
			&A.ptr_[0], &A.idx_[0], &A.val_[0],
			name, opts);
		if(ctx)
			return new direct_solver_A(ctx);
		return 0;
	}

	virtual int solve(const void *b, void *x, size_t nrhs = 1, void *opts = 0) {
		return hj_sparse_direct_solver_A_solve(ctx_, b, x, nrhs, opts);
	}
	~direct_solver_A(void) {
		hj_sparse_direct_solver_A_destroy(ctx_);
	}
protected:
	direct_solver_A(void *ctx):solver_base(ctx){}
};

//! @brief solve AA^Tx=Ab
class direct_solver_AAT : public solver_base
{
public:
	template <typename VAL_TYPE, typename INT_TYPE>
	static direct_solver_AAT *create(
		const csc<VAL_TYPE, INT_TYPE> &A, const char *name = 0, void *opts = 0) {
		return 0;
	}
	virtual int solve(const void *b, void *x, const void *nrhs = 0, void *opts = 0) {
		return 0;
	}
protected:
	direct_solver_AAT(void *ctx):solver_base(ctx){}
};

//! @brief solve Ax=b with iterative method
class iterative_solver_A : public solver_base
{
public:
	template <typename VAL_TYPE, typename INT_TYPE>
	static iterative_solver_A *create(
		const csc<VAL_TYPE, INT_TYPE> &A, const char *name = 0, void *opts = 0) {
		void *ctx = hj_sparse_iterative_solver_A_create(
			sizeof(INT_TYPE), value_type<VAL_TYPE>::size, value_type<VAL_TYPE>::type,
			A.size(1), A.size(2),
			&A.ptr_[0], &A.idx_[0], &A.val_[0],
			name, opts);
		if(ctx)
			return new iterative_solver_A(ctx);
		return 0;
	}
	virtual int solve(const void *b, void *x, size_t nrhs = 1, void *opts = 0) {
		return hj_sparse_iterative_solver_A_solve(ctx_, b, x, nrhs, opts);
	}
	~iterative_solver_A(void) {
		hj_sparse_iterative_solver_A_destroy(ctx_);
	}
protected:
	iterative_solver_A(void *ctx):solver_base(ctx){}
};

}}

#endif
