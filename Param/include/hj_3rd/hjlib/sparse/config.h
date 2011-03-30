#ifndef HJ_SPARSE_CONFIG_H_
#define HJ_SPARSE_CONFIG_H_

#ifdef _WIN32
#  ifdef hj_sparse_solver_EXPORTS
#    define HJ_SPARSE_API __declspec(dllexport)
#  endif
#endif

#ifndef HJ_SPARSE_API
#  define HJ_SPARSE_API
#endif

#endif
