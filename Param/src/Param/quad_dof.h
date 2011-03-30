#ifndef _HJ_QUAD_DOF_H_
#define _HJ_QUAD_DOF_H_

/**
   x: 5x1
   c: 3x2
   E: 1x1
   E_J: 5x1
   E_H: 5x5
 */

void calc_E(double *E, const double *x, const double *c);

void calc_E_J(double *E_J, const double *x, const double *c);

void calc_E_H(double *E_H, const double *x, const double *c);

#endif
