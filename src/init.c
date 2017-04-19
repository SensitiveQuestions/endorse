#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void R2GeoCount (double *dX_lon, double *dX_lat, int *n_X, double *dY_lon,
		 double *dY_lat, int *n_Y, double *Z, int *Store_count);

void R2GeoId (double *dX_lon, double *dX_lat, int *n_X, double *dY_lon,
	      double *dY_lat, int *n_Y, double *Z, int *Store_count);

void R2endorse(int *dY,	int *dT, double *dZ, double *dV, int *village,
	       int *param_indiv, int *n_cat, int *param_village, double *X,
	       double *dS, double *sbeta, double *stau, double *slambda,
	       double *somega2, double *dtheta, double *phi2, double *kappa,
	       double *psi2, double *delta, double *zeta, double *rho2,
	       double *dmu_beta, double *dA0_beta, double *dA0_x, double *dmu_theta,
	       double *dA0_theta, double *dmu_kappa, double *dA0_kappa,
	       double *mu_delta, double *dA0_delta, double *mu_zeta, double *dA0_zeta,
	       double *param_invchisq, int *settings, double *prop, double *betaStore,
	       double *tauStore, double *xStore, double *sStore, double *lambdaStore,
	       double *thetaStore, double *kappaStore, double *deltaStore,
	       double *zetaStore, double *omega2Store, double *phi2Store,
	       double *psi2Store, double *sig2Store, double *rho2Store,
	       double *betaLast, double *tauLast, double *xLast, double *sLast,
	       double *lambdaLast, double *thetaLast, double *kappaLast,
	       double *deltaLast, double *zetaLast, double *omega2Last,
	       double *phi2Last, double *psi2Last, double *sig2Last, double *rho2Last,
	       double *accept_ratio);

static R_NativePrimitiveArgType R2GeoCount_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType R2GeoId_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType R2endorse_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP
};

static const R_CMethodDef CEntries[] = {
  {"R2endorse",  (DL_FUNC) &R2endorse,  64, R2endorse_t},
  {"R2GeoCount", (DL_FUNC) &R2GeoCount,  8, R2GeoCount_t},
  {"R2GeoId",    (DL_FUNC) &R2GeoId,     8, R2GeoId_t},
  {NULL, NULL, 0}
};

void R_init_endorse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
