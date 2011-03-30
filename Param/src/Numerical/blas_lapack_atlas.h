#ifndef HJ_BLAS_LAPACK_ATLAS_H_
#define HJ_BLAS_LAPACK_ATLAS_H_

extern "C" {
#undef small	// <-windows.h <- NpcNdr.h
#include "../hj_3rd/include/lapack/f2c.h"
#undef min
#undef max
#undef abs
#undef dabs
#undef min
#undef max
#undef dmin
#undef dmax
#undef bit_test
#undef bit_clear
#undef bit_set
#undef VOID
#include "../hj_3rd/include/lapack//cblas.h"
#include "../hj_3rd/include/lapack//clapack.h"
}

#endif
