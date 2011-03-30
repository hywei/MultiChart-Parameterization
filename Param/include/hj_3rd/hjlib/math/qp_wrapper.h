#ifndef HJ_QP_WRAPER_H_
#define HJ_QP_WRAPER_H_

#include "conf.h"

#include <zjucad/matrix/matrix.h>

int HJ_MATH_API qp_solve(int nec, const zjucad::matrix::matrix<double> &A, const zjucad::matrix::matrix<double> &B,	// (A^T*X+B)_i: = 0 (when i < nec); > 0 (when i >= nec)
			 const zjucad::matrix::matrix<double> &Aw, const zjucad::matrix::matrix<double> &Bw, const zjucad::matrix::matrix<double> &W,	// min\|diag(W)*(Aw^T*X+Bw)\|^2
			 const zjucad::matrix::matrix<double> &Xmin, const zjucad::matrix::matrix<double> &Xmax,	// bound
			 zjucad::matrix::matrix<double> &X);	// X must be pre-allocated

#endif
