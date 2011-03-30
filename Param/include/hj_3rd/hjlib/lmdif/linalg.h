#ifndef HJ_LINALG_H_
#define HJ_LINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*lmdif_func)(void *ctx, int m, int n, double x[], double fvec[], int *iflag);

void lmdif(lmdif_func func, void *ctx, int m, int n, double x[], double fvec[], double ftol, double xtol, double gtol, int maxfev, double epsfcn,
 double diag[], int mode, double factor, int nprint, int *info, int *nfev, double fjac[],
 int ldfjac,int ipvt[], double qtf[], double wa1[], double wa2[], double wa3[], double wa4[]);

#ifdef __cplusplus
}
#endif

#endif

