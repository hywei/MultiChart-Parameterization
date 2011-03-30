#ifndef HJ_BLAS_LAPACK_H_
#define HJ_BLAS_LAPACK_H_

#ifdef __cplusplus
namespace f2c { extern "C" {
#endif
//#include <f2c.h>
#include "../lapack/f2c.h"
#ifdef __cplusplus
} }
#endif

#undef min
#undef max
#undef abs

#if defined(_WIN32)
#undef small	// <-windows.h <- NpcNdr.h
#endif

#ifdef __cplusplus
namespace clapack {
	using namespace f2c;
#endif
//#include <clapack.h>
	extern "C" {
#include "../lapack/clapack.h"
	}
#ifdef __cplusplus
}
#endif

//#include <cblas.h>
#include "../lapack/cblas.h"

#endif
