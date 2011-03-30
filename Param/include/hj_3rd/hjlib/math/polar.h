#ifndef HJ_POLAR_H_
#define HJ_POLAR_H_

#include "conf.h"

#include <zjucad/matrix/matrix.h>

namespace hj {

class HJ_MATH_API polar3f
{
public:
	polar3f();
	~polar3f();
	// return -1 for bad numerical result
	// others for iteration steps.
	int operator()(zjucad::matrix::matrix<float> &A,	// returned as R
		int reflect = 0, int steps = 5, float eps = 1e-19f) const;
	int operator()(const zjucad::matrix::matrix<float> &A, zjucad::matrix::matrix<float> &R,
		int reflect = 0, int steps = 5, float eps = 1e-19f) const;
private:
	mutable void *ctx_;
};

class HJ_MATH_API polar3d
{
public:
	polar3d();
	~polar3d();
	// return -1 for bad numerical result
	// others for iteration steps.
	int operator()(zjucad::matrix::matrix<double> &A,	// returned as R
		int reflect = 0, int steps = 5, double eps = 1e-19f) const;	// reflect: 0 no, 1 do, 2 auto
	int operator()(const zjucad::matrix::matrix<double> &A, zjucad::matrix::matrix<double> &R,
		int reflect = 0, int steps = 5, double eps = 1e-19f) const;
private:
	mutable void *ctx_;
};

class HJ_MATH_API polar2f
{
public:
	polar2f();
	~polar2f();
	int operator()(zjucad::matrix::matrix<float> &A, int steps = 5, float eps = 1e-19f) const;
private:
	mutable void *ctx_;
};

class HJ_MATH_API polar2d
{
public:
	polar2d();
	~polar2d();
	int operator()(zjucad::matrix::matrix<double> &A, int steps = 5, double eps = 1e-19f) const;
private:
	mutable void *ctx_;
};

}

#endif
