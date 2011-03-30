#ifndef HJ_UMF_PACK_SOLVER_H_
#define HJ_UMF_PACK_SOLVER_H_

#include "sparse.h"

class umf_pack_solver : public hj::sparse::solver
{
public:
	umf_pack_solver();
	bool create(const hj::sparse::spm_csc<double> &A);
	virtual bool set_patten(int rows, const zjucad::matrix::matrix<int> &ptr, const zjucad::matrix::matrix<int> &idx);
	virtual bool set_value (const zjucad::matrix::matrix<double> &val);
	virtual bool solve(const double *b, double *x, int nrhs = 1);
	~umf_pack_solver();
protected:
	void *numeric_, *symbolic_;
	hj::sparse::spm_csc<double> A_;
};

#endif
