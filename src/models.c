#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/*** 
     Bayesian Normal Regression: see Chap.14 of Gelman et al. (2004) 
       both proper and improper priors (and their combinations)
       allowed for beta and sig2. 
       Model: Y_i ~ N(X_i^T beta, sig2)
       Prior: 
         if conjugate = 1,
           Prior for beta: p(beta | sig2) = N(beta | beta0, sig2 * A0)
	   Prior for sig2: p(sig2) = inv-Chi2(sig2 | nu0, s0)
	 if conjugate = 0 (semi-conjugate prior),
	   Prior for beta: p(beta) = N(beta | beta0, A0^{-1})
	   Prior for sig2: p(sig2) = inv-Chi2(sig2 | nu0, s0)
       If conjugate = 1, sig2 is sampled from its marginal
      and beta is sampled from its conditional given sig2.
       If conjugate = 0, beta is updated using its conditional given
      sig2 and sig2 is updated using its conditional given beta.
      In this case, the starting value of sig2 must be provided.
***/
void bNormalReg(double **D,    /* data [X Y] */
		double *beta,  /* coefficients */
		double *sig2,  /* variance */
		int n_samp,    /* sample size */
		int n_cov,     /* # of covariates */
		int addprior,  /* Should prior on beta be incorporated
				  into D? */
		int pbeta,     /* Is prior proper for beta? */
		double *beta0, /* prior mean for normal */
		double **A0,   /* prior precision for normal; can be
				  set to zero to induce improper prior
				  for beta alone
			       */
		int psig2,     /* 0: improper prior for sig2
				  p(sig2|X) \propto 1/sig2
				  1: proper prior for sig2
				  p(sigma2|X) = InvChi2(nu0, s0)
			       */
		double s0,     /* prior scale for InvChi2 */
		int nu0,       /* prior d.f. for InvChi2 */
		int sig2fixed,  /* 1: sig2 fixed, 0: sig2 sampled */
		int conjugate
		) {
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double **mtemp = doubleMatrix(n_cov, n_cov);
  double temp;

  /* slice sampler */
  double f, y;
  double w;
  w = .1;
  int m;
  m = 100;
  
  double L, R;
  int J, K;

  double x;

  /* storage parameters and loop counters */
  int i, j, k;  

  /* if semi-conjugate model, divide y and X by sig2 */
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] /= sqrt(sig2[0]);
	}
      }
    }
  }

  /* read the proper prior for beta as additional data points */
  if (addprior) {
    if (pbeta) {
      dcholdc(A0, n_cov, mtemp);
    } else {
      for (i = 0; i < n_cov; i++)
	for (j = 0; j < n_cov; j++)
	  mtemp[i][j] = 0;
    }
    for(i = 0; i < n_cov; i++) {
      D[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	D[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	D[n_samp+i][j] = mtemp[i][j];
      }
    }
  } 

  /* SS matrix */
  for(j = 0; j <= n_cov; j++) {
    for(k = 0; k <= n_cov; k++) {
      SS[j][k]=0;
    }
  }
  for(i = 0;i < n_samp + n_cov; i++) {
    for(j = 0;j <= n_cov; j++) {
      for(k = 0; k <= n_cov; k++) {
	SS[j][k] += D[i][j]*D[i][k];
      }
    }
  }
  
  /* if semi-conjugate model, y and X are scaled back*/
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] *= sqrt(sig2[0]);
	}
      }
    }
  }
  
  
  /* SWEEP SS matrix */
  for(j = 0; j < n_cov; j++)
    SWP(SS, j, n_cov+1);

  if (pbeta) {
    if (!conjugate) {/* semi-conjugate prior case.
			prior for beta is proper */

      for(j = 0; j < n_cov; j++)
	mean[j] = SS[j][n_cov];

      /* draw beta from its conditional given sig2 */
      for(j = 0; j < n_cov; j++)
	for(k = 0; k < n_cov; k++)
	  V[j][k] = -SS[j][k];




      /* /\*** Debugging ***\/ */
      /* Rprintf("Parameter: mean of coefs\n"); */
      /* PdoubleArray(mean, n_cov); */
      /* Rprintf("Parameter: cov of coefs\n"); */
      /* PdoubleMatrix(V, n_cov, n_cov); */
      /* /\*******\/ */



      rMVN(beta, mean, V, n_cov);

      /* draw sig2 from its conditional given beta */
      /** sum of squared residuals  **/
      SS[n_cov][n_cov] = 0;
      for (i = 0; i < n_samp; i++) {
	temp = 0;
	for (j = 0; j < n_cov; j++) {
	  temp += D[i][j] * beta[j];
	}
	SS[n_cov][n_cov] += (D[i][n_cov] - temp) * (D[i][n_cov] - temp);
      }


      /* /\*** debugging ***\/ */
      /* Rprintf("Sum of squared residuals, %13g\n", SS[n_cov][n_cov]); */
      /* /\*** debugging***\/ */


      /* draw sig2 */
      if (!sig2fixed) {
	if (psig2) {
	  sig2[0] = ( SS[n_cov][n_cov] + nu0 * s0) / rchisq((double)n_samp+nu0);
	} else {
	  if (n_samp > 8) { /* if the inverse Chi-squared distribution is proper,
			       sample sig2 from the truncated inverse Chi-squared */
	    sig2[0] = TruncInvChisq(n_samp - 2, SS[n_cov][n_cov] / (n_samp - 2), 10, 0);

	  } else { /* if the inverse Chi-squared distribution is improper */
	    /* slice sampler code with prior for sig2 being sig2^{-1} with support [0, 10]  */
	    f = pow(sig2[0], - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * sig2[0]));
	    y = runif(0, f);
	    /** "stepping out" procedure **/
	    L = sig2[0] - w * runif(0, 1);
	    if (L < 0) L = 0;
	    R = L + w;
	    if (R > 10) R = 10;
	    J = floor(m * runif(0, 1));
	    K = (m - 1) - J;
	    f = pow(L, - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * L));
	    while (J > 0 && y < f && L > 0) {
	      L = L - w;
	      if (L < 0) {
		L = 0;
		break;
	      }
	      J = J - 1;
	      f = pow(L, - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * L));
	    }
	    f = pow(R, - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * R));
	    while (K > 0 && y < f && R < 10) {
	      R = R + w;
	      if (R > 10) {
		R = 10;
		break;
	      }
	      K = K - 1;
	      f = pow(R, - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * R));
	    }
	    /** "shrinkage" procedure **/
	    do {
	      x = runif(L, R);
	      f = pow(x, - (double)(n_samp / 2)) * exp(- SS[n_cov][n_cov] / (2 * x));
	      if (x < sig2[0]) {
		L = x;
	      } else {
		R = x;
	      }
	    } while (y > f);
	    sig2[0] = x;
	  }

	}
      }
    
    }
  }

  if (!pbeta || conjugate) { /*  conjugate case or improper prior for beta */

    /* draw sig2 from its marginal dist */
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (!sig2fixed) {
      if (psig2) {  /* proper prior for sig2 */
	if (pbeta)   /* proper prior for beta */
	  sig2[0]=(SS[n_cov][n_cov]+nu0*s0)/rchisq((double)n_samp+nu0);
	else        /* improper prior for beta */
	  sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0);
      } else         /* improper prior for sig2 */
	sig2[0]=SS[n_cov][n_cov]/rchisq((double)n_samp-n_cov);
    }

    /* draw beta from its conditional given sig2 */
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2[0];
    rMVN(beta, mean, V, n_cov);

  }

  /* freeing memory */
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/* /\***  Don't use. Obsolete */
/*      Bayesian Normal Regression for "endorse" */
/*        see Chap.14 of Gelman et al. (2004)  */
/*        both proper and improper priors (and their combinations) */
/*        allowed for beta and sig2. */
/* ***\/ */
/* void endorsebNormalReg(double **D,    /\* data [X Y] *\/ */
/* 		       double *beta,  /\* coefficients *\/ */
/* 		       double *sig2,  /\* variance *\/ */
/* 		       int n_samp,    /\* sample size *\/ */
/* 		       int n_cov,     /\* # of covariates *\/ */
/* 		       int addprior,  /\* Should prior on beta be incorporated */
/* 					 into D? *\/ */
/* 		       int pbeta,     /\* Is prior proper for beta? *\/ */
/* 		       double *beta0, /\* prior mean for normal *\/ */
/* 		       double **A0,   /\* prior precision for normal; can be */
/* 					 set to zero to induce improper prior */
/* 					 for beta alone */
/* 				      *\/ */
/* 		       int psig2,     /\* 0: improper prior for sig2 */
/* 					 p(sig2|X) \propto 1/sig2 */
/* 					 1: proper prior for sig2 */
/* 					 p(sigma2|X) = InvChi2(nu0, s0) */
/* 				      *\/ */
/* 		       double s0,     /\* prior scale for InvChi2 *\/ */
/* 		       int nu0,       /\* prior d.f. for InvChi2 *\/ */
/* 		       int sig2fixed  /\* 1: sig2 fixed, 0: sig2 sampled *\/  */
/* 		       ) { */
/*   /\* model parameters *\/ */
/*   double **SS = doubleMatrix(n_cov+1, n_cov+1); /\* matrix folders for SWEEP *\/ */
/*   double *mean = doubleArray(n_cov);            /\* means for beta *\/ */
/*   double **V = doubleMatrix(n_cov, n_cov);      /\* variances for beta *\/ */
/*   double **mtemp = doubleMatrix(n_cov, n_cov); */

/*   /\* storage parameters and loop counters *\/ */
/*   int i, j, k;   */



/*   /\* read the proper prior for beta as additional data points *\/ */
/*   if (addprior) { */
/*     dcholdc(A0, n_cov, mtemp); */
/*     for(i = 0; i < n_cov; i++) { */
/*       D[n_samp+i][n_cov] = 0; */
/*       for(j = 0; j < n_cov; j++) { */
/* 	D[n_samp+i][n_cov] += mtemp[i][j]*beta0[j]; */
/* 	D[n_samp+i][j] = mtemp[i][j]; */
/*       } */
/*     } */
/*   }  */

/* /\*   PdoubleMatrix(D, n_samp + 1, n_cov + 1); *\/ */
  
/*   /\* SS matrix *\/ */
/*   for(j = 0; j <= n_cov; j++) { */
/*     for(k = 0; k <= n_cov; k++) { */
/*       SS[j][k]=0; */
/*     } */
/*   } */
/*   for(i = 0;i < n_samp + n_cov; i++) { */
/*     for(j = 0;j <= n_cov; j++) { */
/*       for(k = 0; k <= n_cov; k++) { */
/* 	SS[j][k] += D[i][j]*D[i][k]; */
/*       } */
/*     } */
/*   } */
  
/* /\*   Rprintf("n_samp %3d \n", n_samp); *\/ */

/* /\*   Rprintf("Sweep matrix %3g \n", SS[n_cov][n_cov]); *\/ */
  
  
/*   /\* SWEEP SS matrix *\/ */
/*   for(j = 0; j < n_cov; j++) */
/*     SWP(SS, j, n_cov+1); */


/* /\*   /\\* check scale parameters *\\/ *\/ */
/* /\*   Rprintf("Prior df %3d, Prior scale %3g, Sweep %3g \n", *\/ */
/* /\* 	  nu0, s0, SS[n_cov][n_cov]); *\/ */


/*   /\* draw sig2 from its marginal dist *\/ */
/*   for(j = 0; j < n_cov; j++) */
/*     mean[j] = SS[j][n_cov]; */
/*   if (!sig2fixed) { */
/*     if (psig2) {  /\* proper prior for sig2 *\/ */
/*       if (pbeta) {  /\* proper prior for beta *\/ */
/* 	sig2[0] = TruncInvChisq((n_samp + nu0), */
/* 				((SS[n_cov][n_cov]+nu0*s0) /  */
/* 				 ((double)(n_samp + nu0))), */
/* 				10, 0); */
/* 	Rprintf("degrees of freedom %3d\n", (n_samp + nu0)); */
/* 	Rprintf("scale %5g\n", ((SS[n_cov][n_cov]+nu0*s0) / ((double)(n_samp + nu0)))); */
/* 	Rprintf("sampled phi2: %5g\n", sig2[0]); */
/* 	/\*** Prior is U[0, 10] ***\/ */
/* /\* 	sig2[0] = TruncInvChisq((n_samp - n_cov), *\/ */
/* /\* 				(SS[n_cov][n_cov] / ((double)(n_samp - n_cov))), *\/ */
/* /\* 				10); *\/ */
/*       } else {       /\* improper prior for beta *\/ */
/* 	sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0); */
/*       } */
/*     } else {        /\* improper prior for sig2 *\/ */
/*       sig2[0]=SS[n_cov][n_cov]/rchisq((double)n_samp-n_cov); */
/*     } */
/*   } */

/*   /\* draw beta from its conditional given sig2 *\/ */
/*   for(j = 0; j < n_cov; j++) */
/*     for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2[0]; */
/*   rMVN(beta, mean, V, n_cov); */

/*   /\* freeing memory *\/ */
/*   free(mean); */
/*   FreeMatrix(SS, n_cov+1); */
/*   FreeMatrix(V, n_cov); */
/*   FreeMatrix(mtemp, n_cov); */
/* } */



/*** 
   A Gibbs Sampler for Binary Probit Regression With and Without
   Marginal Data Augmentation
   
   Marginal Data Augmentation: see p.318 of Imai and van Dyk (2005)
   Journal of Econometrics.
      Prior mean for beta will be set to zero. 
      Improper prior allowed (set A0 to be a matrix of zeros).
***/ 

void bprobitGibbs(int *Y,        /* binary outcome variable */
		  double **X,    /* covariate matrix */
		  double *beta,  /* coefficients */
		  int n_samp,    /* # of obs */ 
		  int n_cov,     /* # of covariates */
		  int prior,     /* Should prior be included in X? */
		  double *beta0, /* prior mean */
		  double **A0,   /* prior precision */
		  int mda,       /* Want to use marginal data augmentation? */ 
		  int n_gen      /* # of gibbs draws */
		  ) {
  
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  
  /* marginal data augmentation */
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	if (!mda)
	  X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /* marginal data augmentation */
    if (mda) sig2 = s0/rchisq((double)nu0);
    
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_cov; j++) 
	dtemp += X[i][j]*beta[j]; 
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+1000,dtemp,1,0);
      X[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    for(i = n_samp;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];

    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);

    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mean, V, n_cov);
 
    /* rescaling the parameters */
    if(mda) {
      for (j = 0; j < n_cov; j++) 
	beta[j] /= sqrt(sig2);
      for (i = 0; i < n_samp; i++)
	X[i][n_cov] /= sqrt(sig2);
    }
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(W);
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/*** 
   A Gibbs Sampler for Ordinal Probit Regression With and Without
   Marginal Data Augmentation
   
   Marginal Data Augmentation for updating coefficients: 
   p.318 of Imai and van Dyk (2005) Journal of Econometrics.
   Prior mean for beta will be set to zero. 
   Improper prior allowed (set A0 to be a matrix of zeros). 
  
   Metropolis-Hasting Blocking Step for updating cutpoints:
   Cowles (1996). Statistics and Computing
   
***/ 

void boprobitMCMC(int *Y,        /* ordinal outcome variable: 0, 1,
				    dots, J-1 */
		  double **X,    /* covariate matrix */
		  double *beta,  /* coefficients */
		  double *tau,   /* J cut points: the first
				    cutpoint is set to 0 and the last
				    cutpoint is set to tau_{J-1}+1000 */
		  int n_samp,    /* # of obs */ 
		  int n_cov,     /* # of covariates */
		  int n_cat,     /* # of categories: J */
		  int prior,     /* Should prior be included in X? */
		  double *beta0, /* prior mean */
		  double **A0,   /* prior precision */
		  int mda,       /* use marginal data augmentation? */
		  int mh,        /* use metropolis-hasting step? */
		  double *prop,  /* J-2 proposal variances for MH step */
		  int *accept,   /* counter for acceptance */
		  int n_gen      /* # of gibbs draws */
		  ) {

  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_samp);           /* means for each obs */
  double *mbeta = doubleArray(n_cov);           /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double *Wmax = doubleArray(n_cat);  /* max of W in each categry: 0, 1,
					 ..., J-1 */
  double *Wmin = doubleArray(n_cat);  /* min of W in each category: 0, 1, 
					 ..., J-1 */
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  double *dvtemp = doubleArray(n_cat); dvtemp[0] = tau[0];
  double **mtemp = doubleMatrix(n_cov, n_cov);
  
  /* marginal data augmentation */
  double sig2; sig2 = 1;
  int nu0; nu0 = 1;
  double s0; s0 = 1;

  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    for (i = 0; i < n_samp; i++){
      mean[i] = 0;
      for (j = 0; j < n_cov; j++) {
	mean[i] += X[i][j]*beta[j]; 
      }
    }

    /* Sampling tau with MH step */
    if (mh) {
      for (j = 1; j < (n_cat-1); j++) {
	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1);
      }
      dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000;
      for (j = 1; j < (n_cat-1); j++) 
	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) -
			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) -
	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) -
	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0));
      for (i = 0; i < n_samp; i++) {
	if (Y[i] == (n_cat-1))  
	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-mean[i], 0, 1, 0, 1) -
	    pnorm(tau[n_cat-2]-mean[i], 0, 1, 0, 1);
	else if (Y[i] > 0) 
	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-mean[i], 0, 1, 1, 0) -
			      pnorm(dvtemp[Y[i]-1]-mean[i], 0, 1, 1, 0)) -
	    log(pnorm(tau[Y[i]]-mean[i], 0, 1, 1, 0) -
		pnorm(tau[Y[i]-1]-mean[i], 0, 1, 1, 0));
      }
      if (unif_rand() < exp(dtemp)) {
	accept[0]++;
	for (j = 1; j < n_cat; j++) {
	  tau[j] = dvtemp[j];
	}
      }
    } 

    /* Sampling the Latent Variable */
    if (!mh) {
      Wmin[0] = tau[0]; Wmax[0] = tau[0]-10;
      for (j = 1; j < n_cat; j++) {
	Wmin[j] = tau[j];
	Wmax[j] = tau[j-1];
      }
    }

    if (mda) /* marginal data augmentation */ 
      sig2 = s0/rchisq((double)nu0);

    for (i = 0; i < n_samp; i++){
      if (Y[i] == 0) {
	W[i] = TruncNorm(mean[i]-1000,0,mean[i],1,0);
      } else {
	  W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean[i],1,0);
      }
      if (!mh) {
	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]);
	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]);
      }
      X[i][n_cov] = W[i]*sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0; i < n_samp+n_cov; i++)
      for(j = 0; j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);
    
    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mbeta[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mbeta, V, n_cov);
    /* rescaling the parameters */
    if (mda) {
      for (j = 0; j < n_cov; j++) 
	beta[j] /= sqrt(sig2);
      for (i = 0; i < n_samp; i++)
	X[i][n_cov] /= sqrt(sig2);
    }

    /* sampling taus without MH-step */
    if (!mh) { 
      for (j = 1; j < n_cat-1; j++) 
	tau[j] = runif(fmax2(tau[j-1], Wmax[j]), 
		       fmin2(tau[j+1], Wmin[j+1]));
      tau[n_cat-1] = tau[n_cat-2] + 1000;
    }
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  
  /* freeing memory */
  FreeMatrix(SS, n_cov+1);
  free(mean);
  free(mbeta);
  FreeMatrix(V, n_cov);
  free(W);
  free(Wmax);
  free(Wmin);
  free(dvtemp);
  FreeMatrix(mtemp, n_cov);
}


/***
     A Gibbs sampler for ordered probit regression for ideal point estimation.
     works only if n_cov = 2, alpha is sampled from the marginal, and beta
     is sampled from the conditional truncated below at 0.
 ***/


void endorseoprobitMCMC(int *Y,        /* ordinal outcome variable: 0, 1,
					   dots, J-1 */
			 double **X,    /* covariate matrix */
			 double *beta,  /* coefficients */
			 double *tau,   /* J cut points: the first
					   cutpoint is set to 0 and the last
					   cutpoint is set to tau_{J-1}+1000 */
			 int n_samp,    /* # of obs */ 
			 int n_cov,     /* # of covariates */
			 int n_cat,     /* # of categories: J */
			 int prior,     /* Should prior be included in X? */
			 double *beta0, /* prior mean */
			 double **A0,   /* prior precision */
			 int mda,       /* use marginal data augmentation? */
			 int mh,        /* use metropolis-hasting step? */
			 double *prop,  /* J-2 proposal variances for MH step */
			 int *accept,   /* counter for acceptance */
			 int n_gen      /* # of gibbs draws */
			 ) {

  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_samp);           /* means for each obs */
  double *mbeta = doubleArray(n_cov);           /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double *Wmax = doubleArray(n_cat);  /* max of W in each categry: 0, 1,
					 ..., J-1 */
  double *Wmin = doubleArray(n_cat);  /* min of W in each category: 0, 1, 
					 ..., J-1 */
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  double *dvtemp = doubleArray(n_cat); dvtemp[0] = tau[0];
  double **mtemp = doubleMatrix(n_cov, n_cov);
  
  /* marginal data augmentation */
  double sig2; sig2 = 1;
  int nu0; nu0 = 1;
  double s0; s0 = 1;

/*   PdoubleMatrix(A0, n_cov, n_cov); */

  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    for (i = 0; i < n_samp; i++){
      mean[i] = 0;
      for (j = 0; j < n_cov; j++) {
	mean[i] += X[i][j]*beta[j]; 
      }
    }

    /* Sampling tau with MH step */
    if (mh) {
      for (j = 1; j < (n_cat-1); j++) {
	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1);
      }
      dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000;
      for (j = 1; j < (n_cat-1); j++) 
	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) -
			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) -
	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) -
	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0));
      for (i = 0; i < n_samp; i++) {
	if (Y[i] == (n_cat-1))  
	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-mean[i], 0, 1, 0, 1) -
	    pnorm(tau[n_cat-2]-mean[i], 0, 1, 0, 1);
	else if (Y[i] > 0) 
	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-mean[i], 0, 1, 1, 0) -
			      pnorm(dvtemp[Y[i]-1]-mean[i], 0, 1, 1, 0)) -
	    log(pnorm(tau[Y[i]]-mean[i], 0, 1, 1, 0) -
		pnorm(tau[Y[i]-1]-mean[i], 0, 1, 1, 0));
      }
      if (unif_rand() < exp(dtemp)) {
	accept[0]++;
	for (j = 1; j < n_cat; j++) {
	  tau[j] = dvtemp[j];
	}
      }
    } 

    /* Sampling the Latent Variable */
    if (!mh) {
      Wmin[0] = tau[0]; Wmax[0] = tau[0]-10;
      for (j = 1; j < n_cat; j++) {
	Wmin[j] = tau[j];
	Wmax[j] = tau[j-1];
      }
    }

    if (mda) /* marginal data augmentation */ 
      sig2 = s0/rchisq((double)nu0);

    for (i = 0; i < n_samp; i++){
      if (Y[i] == 0) {
	W[i] = TruncNorm(mean[i]-1000,0,mean[i],1,0);
      } else if (Y[i] < 0) {
	W[i] = mean[i] + norm_rand();
      } else {
	W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean[i],1,0);
      }
      if (!mh) {
	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]);
	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]);
      }
      X[i][n_cov] = W[i]*sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;

    for(i = 0; i < n_samp+n_cov; i++)
      for(j = 0; j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];


/*     PdoubleMatrix(X, n_samp + n_cov, n_cov); */
/*     Rprintf("Oprobit Sweep matrix %3g \n", SS[n_cov][n_cov]); */

    
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);
    
/*     /\* check scale parameters *\/ */
/*     Rprintf("Prior df %3g, Prior scale %3g, Sweep %3g \n", */
/* 	    nu0, s0, SS[n_cov][n_cov]); */


    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mbeta[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++)
	V[j][k]=-SS[j][k]*sig2;


/*     Rprintf("   Variance-covariance matrix\n"); */
/*     for (j = 0; j < n_cov; j++) { */
/*       for (k = 0; k < n_cov; k++) */
/* 	Rprintf("   %5g", V[j][k]); */
/*       Rprintf("\n"); */
/*     } */

/*     Rprintf("Mean of alpha                  %14g\n", mbeta[0]); */
/*     Rprintf("Variance of the marginal alpha %14g\n", V[0][0]); */

    beta[0] = mbeta[0] + norm_rand()*sqrt(V[0][0]);
/*     Rprintf("   Sampled alpha: %5g\n", beta[0]); */

/*     Rprintf("Sampled alpha                  %14g\n", beta[0]); */


/*     Rprintf("Mean of beta                   %14g\n", (mbeta[1] + (V[0][1] / V[0][0])*(beta[0] - mbeta[0]))); */
/*     Rprintf("Variance of conditional beta   %14g\n", (V[1][1] - (V[0][1]*V[0][1] / V[0][0]))); */
    
    beta[1] = TruncNorm(0, 10000,
			(mbeta[1] + (V[0][1] / V[0][0])*(beta[0] - mbeta[0])),
			V[1][1] - (V[0][1]*V[0][1] / V[0][0]), 0);

/*     Rprintf("Sampled beta                   %14g\n", beta[1]); */

/*     Rprintf("sig2                           %14g\n", sig2); */

/*     Rprintf("   Sampled beta: %5g\n", beta[1]); */

/*     rMVN(beta, mbeta, V, n_cov); */

    /* rescaling the parameters */
    if (mda) {
      for (j = 0; j < n_cov; j++) 
	beta[j] /= sqrt(sig2);
      for (i = 0; i < n_samp; i++)
	X[i][n_cov] /= sqrt(sig2);
    }

    /* sampling taus without MH-step */
    if (!mh) { 
      for (j = 1; j < n_cat-1; j++) {
	tau[j] = runif(fmax2(tau[j-1], Wmax[j]), 
		       fmin2(tau[j+1], Wmin[j+1]));
      }
      tau[n_cat-1] = tau[n_cat-2] + 1000;
    }
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  
  /* freeing memory */
  FreeMatrix(SS, n_cov+1);
  free(mean);
  free(mbeta);
  FreeMatrix(V, n_cov);
  free(W);
  free(Wmax);
  free(Wmin);
  free(dvtemp);
  FreeMatrix(mtemp, n_cov);
}



/*** 
   A Random Walk Metroplis Sampler for Binomial and Multinomial
   Logistic Regression with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
***/
void logitMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		double **X,    /* (N x K) covariate matrix */
		double *beta,  /* (K(J-1)) stacked coefficient vector */
		int n_samp,    /* # of obs */
		int n_dim,     /* # of categories, J-1 */
		int n_cov,     /* # of covariates, K */
		double *beta0, /* (K(J-1)) prior mean vector */
		double **A0,   /* (K(J-1) x K(J-1)) prior precision */
		double *Var,   /* K(J-1) proposal variances */
		int n_gen,     /* # of MCMC draws */
		int *counter   /* # of acceptance for each parameter */
		) {
  
  int i, j, k, main_loop;
  double numer, denom;
  double *sumall = doubleArray(n_samp); 
  double *sumall1 = doubleArray(n_samp);
  double *prop = doubleArray(n_dim*n_cov);
  double **Xbeta = doubleMatrix(n_samp, n_dim);
  double **Xbeta1 = doubleMatrix(n_samp, n_dim);

  for (j = 0; j < n_cov*n_dim; j++)
    prop[j] = beta[j];
  for (i = 0; i < n_samp; i++) {
    sumall[i] = 1.0; 
    for (j = 0; j < n_dim; j++) {
      Xbeta[i][j] = 0;
      for (k = 0; k < n_cov; k++) 
	Xbeta[i][j] += X[i][k]*beta[j*n_cov+k];
      Xbeta1[i][j] = Xbeta[i][j];
      sumall[i] += exp(Xbeta[i][j]);
    }
    sumall1[i] = sumall[i];
  }

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    for (j = 0; j < n_dim; j++)
      for (k = 0; k < n_cov; k++) {
	/** Sample from the proposal distribution **/
	prop[j*n_cov+k] = beta[j*n_cov+k] + 
	  norm_rand()*sqrt(Var[j*n_cov+k]);
      
      /** Calculating the ratio (log scale) **/
      /* prior */
      numer = dMVN(prop, beta0, A0, n_cov*n_dim, 1);
      denom = dMVN(beta, beta0, A0, n_cov*n_dim, 1);   
      /* likelihood */
      for (i = 0; i < n_samp; i++) {
	Xbeta1[i][j] = Xbeta[i][j] - X[i][k]*(beta[j*n_cov+k]-prop[j*n_cov+k]);
	if (Y[i] > 0) {
	  denom += Xbeta[i][Y[i]-1];
	  numer += Xbeta1[i][Y[i]-1];
	} 
	sumall1[i] += (exp(Xbeta1[i][j])-exp(Xbeta[i][j]));
	numer -= log(sumall1[i]);
	denom -= log(sumall[i]);
      }
      
      /** Rejection **/
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	counter[j*n_cov+k]++;
	beta[j*n_cov+k] = prop[j*n_cov+k];
	for (i = 0; i < n_samp; i++) {
	  sumall[i] = sumall1[i];
	  Xbeta[i][j] = Xbeta1[i][j];
	}
      }
    }
  }
  
  free(prop);
  free(sumall);
  free(sumall1);
  FreeMatrix(Xbeta, n_samp);
  FreeMatrix(Xbeta1, n_samp);
} /* end of logitMetro */



/* /\*** */
/*    A Standard Gibbs Sampler for Normal Mixed Effects Regression */

/*    MODEL: Y_i = X_i \beta + Z_i \gamma_i + \epsilon_i */
/*           where */
/*           \epsilon_i  \sim N(0, \sigma^2 I_{n_i}) */
/*           \gamma_i \sim N(0, \Psi^{-1}) */
/*           and i indexes groups. */
/*    PRIOR: p(\beta|X, Z) = N(beta_0, A_0^{-1}) */
/*           p(\Psi^{-1}|X, Z) = Wish(\tau_0, T_0) */
/*           p(\sigma^2|X, Z) = InvChi2(\nu_0, s_0) */

/* 	  standard diffuse improper priro (imp = 1 && A0 = 0) */
/* 	  p(\beta, log\sigma, \Psi|X, Z) \propto |\Psi|^{-1/2} */
/* 	  The posterior is proper so long as (see van Dyk, 2000, JCGS) */
/*              n_samp > n_fixed + n_random */
/*              n_grp > 2 n_random - 1 */
/* ***\/ */

/* void bNormalMixedGibbs(double *Y,       /\* outcome variable *\/ */
/* 		       double **X,      /\* model matrix for fixed */
/* 					   effects *\/ */
/* 		       double ***Zgrp,  /\* model matrix for random */
/* 					   effects organized by */
/* 					   grous *\/ */
/* 		       int *grp,        /\* group indicator: 0, 1, 2,... *\/ */
/* 		       double *beta,    /\* fixed effects coefficients *\/ */
/* 		       double **gamma,  /\* random effects coefficients *\/ */
/* 		       double *sig2,    /\* variance parameter *\/ */
/* 		       double **Psi,    /\* precision matrix for random */
/* 					   effects *\/ */
/* 		       int n_samp,      /\* # of obs *\/ */
/* 		       int n_fixed,     /\* # of fixed effects *\/ */
/* 		       int n_random,    /\* # of random effects *\/ */
/* 		       int n_grp,       /\* # of groups *\/ */
/* 		       int prior,       /\* include prior for fixed effects in X? *\/ */
/* 		       double *beta0,   /\* prior mean *\/ */
/* 		       double **A0,     /\* prior precision *\/ */
/* 		       int imp,         /\* use standard improper prior */
/* 					   for sig2 and Psi */
/* 					   (see above; beta can be proper) *\/ */
/* 		       int nu0,         /\* prior df for sig2 *\/ */
/* 		       double s0,       /\* prior scale for sig2 *\/ */
/* 		       int tau0,        /\* prior df for Psi *\/ */
/* 		       double **T0,     /\* prior scale for Psi *\/ */
/* 		       int n_gen        /\* # of gibbs draws *\/ */
/* 		       ) { */
  
/*   double *gamma0 = doubleArray(n_random);           /\* prior mean for gamma *\/ */
/*   double **V = doubleMatrix(n_fixed, n_fixed);      /\* variances for beta *\/ */
/*   double **mtemp = doubleMatrix(n_random, n_random); */
/*   double **mtemp1 = doubleMatrix(n_random, n_random); */

/*   /\* storage parameters and loop counters *\/ */
/*   int i, j, k, l, main_loop; */
/*   int *vitemp = intArray(n_grp); */
  
/*   /\* read the prior as additional data points *\/ */
/*   if (prior) { */
/*     dcholdc(A0, n_fixed, V); */
/*     for(i = 0; i < n_fixed; i++) { */
/*       X[n_samp+i][n_fixed] = 0; */
/*       for(j = 0; j < n_fixed; j++) { */
/* 	X[n_samp+i][n_fixed] += V[i][j]*beta0[j]; */
/* 	X[n_samp+i][j] = V[i][j]; */
/*       } */
/*     } */
/*   } */

/*   for (j = 0; j < n_random; j++) */
/*     gamma0[j] = 0; */

/*   /\* Gibbs Sampler! *\/ */
/*   for(main_loop = 1; main_loop <= n_gen; main_loop++){ */
/*     /\** STEP 1: Sample Fixed Effects Given Random Effects */
/*                 Also Sample Variance Parameter **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++) { */
/*       X[i][n_fixed] = Y[i]; */
/*       for (j = 0; j < n_random; j++) */
/* 	X[i][n_fixed] -= Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j]; */
/*       vitemp[grp[i]]++; */
/*     } */
/*     if (imp) */
/*       bNormalReg(X, beta, sig2, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1, */
/* 		 1, 0); */
/*     else */
/*       bNormalReg(X, beta, sig2, n_samp, n_fixed, 0, 1, beta0, A0, 1, s0, */
/* 		 nu0, 0); */

/*     /\** STEP 2: Update Random Effects Given Fixed Effects **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++) { */
/*       Zgrp[grp[i]][vitemp[grp[i]]][n_random] = Y[i]; */
/*       for (j = 0; j < n_fixed; j++) */
/* 	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j]; */
/*       vitemp[grp[i]]++; */
/*     } */
/*     for (j = 0; j < n_grp; j++) */
/*       bNormalReg(Zgrp[j], gamma[j], sig2, vitemp[j], n_random, */
/* 		 1, 1, gamma0, Psi, 0, 0, 1, 1); */

/*     /\** STEP 3: Update Covariance Matrix Given Random Effects **\/ */
/*     for (j = 0; j < n_random; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	if (imp) */
/* 	  mtemp[j][k] = 0; */
/* 	else */
/* 	  mtemp[j][k] = T0[j][k]; */
/*     for (j = 0; j < n_grp; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	for (l = 0; l < n_random; l++) */
/* 	  mtemp[k][l] += gamma[j][k]*gamma[j][l]; */
/*     dinv(mtemp, n_random, mtemp1); */
/*     if (imp) */
/*       rWish(Psi, mtemp1, n_grp-n_random-1, n_random); */
/*     else */
/*       rWish(Psi, mtemp1, tau0+n_grp, n_random); */

/*     R_CheckUserInterrupt(); */
/*   } /\* end of Gibbs sampler *\/ */

/*   /\* freeing memory *\/ */
/*   free(gamma0); */
/*   FreeMatrix(V, n_fixed); */
/*   FreeMatrix(mtemp, n_random); */
/*   FreeMatrix(mtemp1, n_random); */
/*   free(vitemp); */
/* } */



/* /\***  */
/*    A Standard Gibbs Sampler for Binary Probit Mixed Effects Regression */

/*    MODEL: Y_{ij} = 1 if W_{ij} > 0  */
/*                  = 0 otherwise. */
/*           W_i = X_i \beta + Z_i \gamma_i + \epsilon_i  */
/*           where  */
/*           \epsilon_i  \sim N(0, I_{n_i}) */
/*           \gamma_i \sim N(0, \Psi^{-1}) */
/*           and i indexes groups. */
/*    PRIOR: p(\beta|X, Z) = N(beta_0, A_0^{-1}) */
/*           p(\Psi^{-1}|X,Z) = Wish(\tau_0, T_0) */
/*    see the docs for bprobitGibbs for the implmentation of marginal */
/*           data augmentation for fixed effects coefficients        */
/* ***\/  */

/* void bprobitMixedGibbs(int *Y,          /\* binary outcome variable *\/ */
/* 		       double **X,      /\* model matrix for fixed */
/* 					   effects *\/ */
/* 		       double ***Zgrp,  /\* model matrix for random */
/* 					   effects organized by grous *\/ */
/* 		       int *grp,        /\* group indicator: 0, 1, 2,... *\/ */
/* 		       double *beta,    /\* fixed effects coefficients *\/ */
/* 		       double **gamma,  /\* random effects coefficients *\/ */
/* 		       double **Psi,    /\* precision matrix for random */
/* 					   effects *\/ */
/* 		       int n_samp,      /\* # of obs *\/  */
/* 		       int n_fixed,     /\* # of fixed effects *\/ */
/* 		       int n_random,    /\* # of random effects *\/ */
/* 		       int n_grp,       /\* # of groups *\/ */
/* 		       int prior,       /\* include prior for fixed effects in X? *\/ */
/* 		       double *beta0,   /\* prior mean *\/ */
/* 		       double **A0,     /\* prior precision *\/ */
/* 		       int tau0,        /\* prior df *\/ */
/* 		       double **T0,     /\* prior scale *\/ */
/* 		       int n_gen        /\* # of gibbs draws *\/ */
/* 		       ) { */
  
/*   double *gamma0 = doubleArray(n_random);           /\* prior mean for gamma *\/ */
/*   double **V = doubleMatrix(n_fixed, n_fixed);      /\* variances for beta *\/ */
/*   double *W = doubleArray(n_samp); */
/*   double **mtemp = doubleMatrix(n_random, n_random); */
/*   double **mtemp1 = doubleMatrix(n_random, n_random); */

/*   /\* storage parameters and loop counters *\/ */
/*   int i, j, k, l, main_loop;   */
/*   int *vitemp = intArray(n_grp); */
/*   double dtemp0, dtemp1; */
/*   double *vdtemp = doubleArray(1); */
/*   vdtemp[0] = 1.0; */
  
/*   /\* read the prior as additional data points *\/ */
/*   if (prior) { */
/*     dcholdc(A0, n_fixed, V); */
/*     for(i = 0; i < n_fixed; i++) { */
/*       X[n_samp+i][n_fixed] = 0; */
/*       for(j = 0; j < n_fixed; j++) { */
/* 	X[n_samp+i][n_fixed] += V[i][j]*beta0[j]; */
/* 	X[n_samp+i][j] = V[i][j]; */
/*       } */
/*     } */
/*   } */

/*   for (j = 0; j < n_random; j++) */
/*     gamma0[j] = 0; */

/*   /\* Gibbs Sampler! *\/ */
/*   for(main_loop = 1; main_loop <= n_gen; main_loop++){ */
/*     /\** STEP 1: Sample Latent Variable **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++){ */
/*       dtemp0 = 0; dtemp1 = 0; */
/*       for (j = 0; j < n_fixed; j++)  */
/* 	dtemp0 += X[i][j]*beta[j];  */
/*       for (j = 0; j < n_random; j++) */
/* 	dtemp1 += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j]; */
/*       if(Y[i] == 0)  */
/* 	W[i] = TruncNorm(dtemp0+dtemp1-1000,0,dtemp0+dtemp1,1,0); */
/*       else  */
/* 	W[i] = TruncNorm(0,dtemp0+dtemp1+1000,dtemp0+dtemp1,1,0); */
/*       X[i][n_fixed] = W[i]-dtemp1; */
/*       vitemp[grp[i]]++; */
/*     } */
/*     /\** STEP 2: Sample Fixed Effects Given Random Effects **\/ */
/*     bNormalReg(X, beta, vdtemp, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1, */
/* 	       1, 1); */

/*     /\** STEP 3: Update Random Effects Given Fixed Effects **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++) { */
/*       Zgrp[grp[i]][vitemp[grp[i]]][n_random] = W[i]; */
/*       for (j = 0; j < n_fixed; j++)  */
/* 	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j];  */
/*       vitemp[grp[i]]++; */
/*     } */
/*     for (j = 0; j < n_grp; j++) */
/*       bNormalReg(Zgrp[j], gamma[j], vdtemp, vitemp[j], n_random, */
/* 		 1, 1, gamma0, Psi, 0, 0, 1, 1); */

/*     /\** STEP 4: Update Covariance Matrix Given Random Effects **\/ */
/*     for (j = 0; j < n_random; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	mtemp[j][k] = T0[j][k]; */
/*     for (j = 0; j < n_grp; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	for (l = 0; l < n_random; l++) */
/* 	  mtemp[k][l] += gamma[j][k]*gamma[j][l]; */
/*     dinv(mtemp, n_random, mtemp1); */
/*     rWish(Psi, mtemp1, tau0+n_grp, n_random); */

/*     R_CheckUserInterrupt(); */
/*   } /\* end of Gibbs sampler *\/ */

/*   /\* freeing memory *\/ */
/*   free(W); */
/*   free(vdtemp); */
/*   free(vitemp); */
/*   free(gamma0); */
/*   FreeMatrix(V, n_fixed); */
/*   FreeMatrix(mtemp, n_random); */
/*   FreeMatrix(mtemp1, n_random); */
/* } /\* end of mixed effects probit *\/ */

/**
   A Random Walk Metroplis Sampler for Binomial Logistic Mixed Effects 
   Regression with Independent Normal Prior and Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 

**/

void BinomLogitMixed(int *Y,          /* outcome variable: 0, 1, ..., J */
		     double **X,      /* (N x K) covariate matrix for
					 fixed effects */
		     double **Z,      /* (N x L) covariate matrix for 
					 random effects */
		     int *grp,        /* group indicator, 0, 1, ..., G-1 */
		     double *beta,    /* K coefficients for fixed effects */
		     double **gamma,  /* (G x L) matrix of random effects */
		     double **Psi,    /* LxL precision matrix for random effecs */
		     int n_samp,      /* # of obs */
		     int J,           /* size of binomial, J */
		     int n_fixed,     /* # of fixed effects, K */
		     int n_random,    /* # of random effects, L */
		     int n_grp,       /* # of groups, G */
		     double *beta0,   /* K dimensional prior mean vector */
		     double **A0,     /* (K x K) prior precision */
		     int tau0,        /* prior df for Psi */
		     double **T0,     /* prior scale for Psi */
		     double *tune_fixed,  /* K proposal variances */
		     double *tune_random, /* tuning constant for random effects of each group */
		     int n_gen,        /* # of MCMC draws */
		     int *acc_fixed,   /* # of acceptance for fixed effects */
		     int *acc_random   /* # of acceptance for random effects */
		     ) {
  
  int i, j, k, main_loop;
  double numer, denom;
  /* proposal values */
  double *beta1 = doubleArray(n_fixed);
  double *gamma1 = doubleArray(n_random);
  /* prior for gamma = 0 */
  double *gamma0 = doubleArray(n_random);
  /* data holders */
  double *Xbeta = doubleArray(n_samp);
  double *Xbeta1 = doubleArray(n_samp);
  double *Zgamma = doubleArray(n_samp);
  double *Zgamma1 = doubleArray(n_samp);
  /* matrix holders */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed; j++)
    beta1[j] = beta[j];

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /** initializing Xbeta and Zgamma **/
  for (i = 0; i < n_samp; i++) {
    Xbeta[i] = 0; Zgamma[i] = 0;
    for (j = 0; j < n_fixed; j++) { 
      Xbeta[i] += X[i][j] * beta[j];
    }
    Xbeta1[i] = Xbeta[i];
    for (j = 0; j < n_random; j++) {
      Zgamma[i] += Z[i][j] * gamma[grp[i]][j];
    }
    Zgamma1[i] = Zgamma[i];
  }

  /** MCMC Sampler starts here **/
  for (main_loop = 0; main_loop < n_gen; main_loop++) {

    /** STEP 1: Update Each Fixed Effect **/
    for (j = 0; j < n_fixed; j++) {
      /* Sample from the proposal distribution */
      beta1[j] = beta[j] + norm_rand() * sqrt(tune_fixed[j]);
      /* Calculating the ratio (log scale) */
      /* prior */
      numer = dMVN(beta1, beta0, A0, n_fixed, 1);
      denom = dMVN(beta, beta0, A0, n_fixed, 1);   
      /* likelihood */
      for (i = 0; i < n_samp; i++) {
	Xbeta1[i] = Xbeta[i] - X[i][j] * (beta[j] - beta1[j]);
	denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
	numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta1[i]-Zgamma[i])), 1);
      }
      /* Rejection */
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	acc_fixed[j]++;
	beta[j] = beta1[j];
	for (i = 0; i < n_samp; i++) {
	  Xbeta[i] = Xbeta1[i];
	}
      }
    }
 
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    dinv(Psi, n_random, mtemp);
    for (i = 0; i < n_random; i++)
      for (j = 0; j < n_random; j++)
	mtemp[i][j] *= tune_random[j];
    for (j = 0; j < n_grp; j++) {
      rMVN(gamma1, gamma[j], mtemp, n_random);
      /* Calculating the ratio (log scale) */
      /* prior */
      numer = dMVN(gamma1, gamma0, Psi, n_random, 1);
      denom = dMVN(gamma[j], gamma0, Psi, n_random, 1); 
      /* likelihood for group j */
      for (i = 0; i < n_samp; i++) {
	if (grp[i] == j) {
	  for (k = 0; k < n_random; k++)
	    Zgamma1[i] = Zgamma[i] - Z[i][k]*(gamma[j][k]-gamma1[k]);
	  denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
	  numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma1[i])), 1);
	}
      }
      /* Rejection */
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	acc_random[j]++;
	for (k = 0; k < n_random; k++)
	  gamma[j][k] = gamma1[k];
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == j) {
	    Zgamma[i] = Zgamma1[i];
	  }      
	}
      }
    }
    
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (i = 0; i < n_grp; i++)
      for (j = 0; j < n_random; j++)
	for (k = 0; k < n_random; k++)
	  mtemp[j][k] += gamma[i][j] * gamma[i][k];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, tau0+n_grp, n_random);
  }

  /* freeing memory */
  free(beta1);
  free(gamma1);
  free(gamma0);
  free(Xbeta);
  free(Xbeta1);
  free(Zgamma);
  free(Zgamma1);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */


/*** 
   A Random Walk Metroplis Sampler for Binomial and Multinomial
   Logistic Mixed Effects Regression with Independent Normal Prior and
   Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 
***/
void logitMixedMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		     double **X,    /* (N x K) covariate matrix for
				       fixed effects */
		     double ***Z,     /* covariates for random effects
					 organized by groups */
		     int *grp,        /* group indicator, 0, 1, ...,
					 G-1 */
		     double *beta,    /* (K(J-1)) stacked coefficient
					 vector for fixed effects */
		     double ***gamma, /*(J-1)G L array of random
					 effects organized by
					 equations and then by
					 groups */
		     double ***Psi,   /* LxL precision matrix for
					 random effecs for each equation */
		     int n_samp,      /* # of obs */
		     int n_dim,       /* # of categories, J-1 */
		     int n_fixed,     /* # of fixed effects, K */
		     int n_random,    /* # of random effects, L */
		     int n_grp,       /* # of groups, G */
		     double *beta0,   /* (K(J-1)) prior mean vector */
		     double **A0,     /* (K(J-1) x K(J-1)) prior
				       precision */
		     int tau0,        /* prior df for Psi */
		     double **T0,     /* prior scale for Psi */
		     double *tune_fixed,   /* K(J-1) proposal variances */
		     double *tune_random,  /* tuning constant for random
					      effects of each random effect */
		     int n_gen,        /* # of MCMC draws */
		     int *acc_fixed,   /* # of acceptance for fixed effects */
		     int *acc_random   /* # of acceptance for random effects */
		     ) {
  
  int i, j, k, l, main_loop;
  int *vitemp = intArray(n_grp);
  double numer, denom;
  double *sumall = doubleArray(n_samp); 
  double *sumall1 = doubleArray(n_samp);
  double *propb = doubleArray(n_dim*n_fixed);
  double *propg = doubleArray(n_random);
  double *gamma0 = doubleArray(n_random);
  double **Xbeta = doubleMatrix(n_samp, n_dim);
  double **Xbeta1 = doubleMatrix(n_samp, n_dim);
  double **Zgamma = doubleMatrix(n_samp, n_dim);
  double **Zgamma1 = doubleMatrix(n_samp, n_dim);
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed*n_dim; j++)
    propb[j] = beta[j];
  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;
  for (j = 0 ; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    sumall[i] = 1.0; 
    for (j = 0; j < n_dim; j++) {
      Xbeta[i][j] = 0; Zgamma[i][j] = 0;
      for (k = 0; k < n_fixed; k++) 
	Xbeta[i][j] += X[i][k]*beta[j*n_fixed+k];
      Xbeta1[i][j] = Xbeta[i][j];
      for (k = 0; k < n_random; k++)
	Zgamma[i][j] += Z[j][vitemp[grp[i]]][k]*gamma[j][grp[i]][k];
      sumall[i] += exp(Xbeta[i][j] + Zgamma[i][j]);
      Zgamma1[i][j] = Zgamma[i][j];
    }
    sumall1[i] = sumall[i];
    vitemp[grp[i]]++;
  }

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
     /** STEP 1: Update Fixed Effects Given Random Effects **/
    for (j = 0; j < n_dim; j++)
      for (k = 0; k < n_fixed; k++) {
	/** Sample from the proposal distribution **/
	propb[j*n_fixed+k] = beta[j*n_fixed+k] + 
	  norm_rand()*sqrt(tune_fixed[j*n_fixed+k]);
	/** Calculating the ratio (log scale) **/
	/* prior */
	numer = dMVN(propb, beta0, A0, n_fixed*n_dim, 1);
	denom = dMVN(beta, beta0, A0, n_fixed*n_dim, 1);   
	/* likelihood */
	for (i = 0; i < n_samp; i++) {
	  Xbeta1[i][j] = Xbeta[i][j] - X[i][k]*(beta[j*n_fixed+k]-propb[j*n_fixed+k]);
	  if (Y[i] > 0) {
	    denom += (Xbeta[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	    numer += (Xbeta1[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	  } 
	  sumall1[i] += (exp(Xbeta1[i][j] + Zgamma[i][Y[i]-1]) - 
			 exp(Xbeta[i][j] + Zgamma[i][Y[i]-1]));
	  numer -= log(sumall1[i]);
	  denom -= log(sumall[i]);
	}
	/** Rejection **/
	if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	  acc_fixed[j*n_fixed+k]++;
	  beta[j*n_fixed+k] = propb[j*n_fixed+k];
	  for (i = 0; i < n_samp; i++) {
	    sumall[i] = sumall1[i];
	    Xbeta[i][j] = Xbeta1[i][j];
	  }
	}
      }
 
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_dim; j++) {
      dinv(Psi[j], n_random, mtemp);
      for (i = 0; i < n_random; i++)
	for (k = 0; k < n_random; k++)
	  mtemp[i][k] *= tune_random[j];
      for (k = 0; k < n_grp; k++) {
	rMVN(propg, gamma[j][k], mtemp, n_random);
	/** Calculating the ratio (log scale) **/
	/* prior */
	numer = dMVN(propg, gamma0, Psi[j], n_random, 1);
	denom = dMVN(gamma[j][k], gamma0, Psi[j], n_random, 1); 
 	/* likelihood */
	for (l = 0; l < n_grp; l++)
	  vitemp[l] = 0;
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == k)
	    for (l = 0; l < n_random; l++)
	      Zgamma1[i][j] -= Z[k][vitemp[k]][l]*(gamma[j][k][l]-propg[l]);
	  vitemp[grp[i]]++;
	  if (Y[i] > 0) {
	    denom += (Xbeta[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	    numer += (Xbeta[i][Y[i]-1] + Zgamma1[i][Y[i]-1]);
	  } 
	  sumall1[i] += (exp(Xbeta[i][j] + Zgamma1[i][j]) -
			 exp(Xbeta[i][j] + Zgamma[i][j]));
	  numer -= log(sumall1[i]);
	  denom -= log(sumall[i]);
	}
	/* Rejection */
	if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	  acc_random[j*n_grp+k]++;
	  for (l = 0; l < n_random; l++)
	    gamma[j][k][l] = propb[l];
	  for (i = 0; i < n_samp; i++) {
	    sumall[i] = sumall1[i];
	    Zgamma[i][j] = Zgamma1[i][j];
	  }      
	}
      }
    }
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_dim; j++) {
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] = T0[k][l];
      for (i = 0; i < n_grp; i++)
	for (k = 0; k < n_random; k++)
	  for (l = 0; l < n_random; l++)
	    mtemp[k][l] += gamma[j][i][k]*gamma[j][i][l];
      dinv(mtemp, n_random, mtemp1);
      rWish(Psi[j], mtemp1, tau0+n_grp, n_random);
    }
  }

  /* freeing memory */
  free(vitemp);
  free(sumall);
  free(sumall1);
  free(propb);
  free(propg);
  free(gamma0);
  FreeMatrix(Xbeta, n_samp);
  FreeMatrix(Xbeta1, n_samp);
  FreeMatrix(Zgamma, n_samp);
  FreeMatrix(Zgamma1, n_samp);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */




/* /\***  */
/*    A Standard Gibbs Sampler for Ordinal Probit Mixed Effects Regression */
/* ***\/  */

/* void boprobitMixedMCMC(int *Y,          /\* binary outcome variable *\/ */
/* 		       double **X,      /\* model matrix for fixed */
/* 					   effects: the last column */
/* 					   contains the starting values for W-Zgamma *\/ */
/* 		       double ***Zgrp,  /\* model matrix for random */
/* 					   effects organized by grous *\/ */
/* 		       int *grp,        /\* group indicator: 0, 1, 2,... *\/ */
/* 		       double *beta,    /\* fixed effects coefficients *\/ */
/* 		       double **gamma,  /\* random effects coefficients *\/ */
/* 		       double *tau,     /\* cutpoints *\/ */
/* 		       double **Psi,    /\* precision matrix for random */
/* 					   effects *\/ */
/* 		       int n_samp,      /\* # of obs *\/  */
/* 		       int n_cat,       /\* number of categories *\/ */
/* 		       int n_fixed,     /\* # of fixed effects *\/ */
/* 		       int n_random,    /\* # of random effects *\/ */
/* 		       int n_grp,       /\* # of groups *\/ */
/* 		       int prior,       /\* include prior for fixed effects in X? *\/ */
/* 		       double *beta0,   /\* prior mean *\/ */
/* 		       double **A0,     /\* prior precision *\/ */
/* 		       int tau0,        /\* prior df *\/ */
/* 		       double **T0,     /\* prior scale *\/ */
/* 		       int mh,          /\* metropolis-hastings step */
/* 					   for cutpoints? *\/ */
/* 		       double *prop,    /\* proposal variance for MH */
/* 					   step *\/ */
/* 		       int *accept,     /\* counter for acceptance *\/ */
/* 		       int n_gen        /\* # of gibbs draws *\/ */
/* 		       ) { */
  
/*   double *gamma0 = doubleArray(n_random);         /\* prior mean for gamma *\/ */
/*   double *Xbeta = doubleArray(n_samp);            /\* X beta *\/ */
/*   double *Zgamma = doubleArray(n_samp); */
/*   double **V = doubleMatrix(n_fixed, n_fixed);    /\* variances for beta *\/ */
/*   double *W = doubleArray(n_samp); */
/*   double *Wmax = doubleArray(n_cat);  /\* max of W in each categry: 0, 1, */
/* 					 ..., J-1 *\/ */
/*   double *Wmin = doubleArray(n_cat);  /\* min of W in each category: 0, 1,  */
/* 					 ..., J-1 *\/ */
/*   double **mtemp = doubleMatrix(n_random, n_random); */
/*   double **mtemp1 = doubleMatrix(n_random, n_random); */

/*   /\* storage parameters and loop counters *\/ */
/*   int i, j, k, l, main_loop;   */
/*   int *vitemp = intArray(n_grp); */
/*   double dtemp; */
/*   double *vdtemp = doubleArray(1); */
/*   double *dvtemp = doubleArray(n_cat); */
/*   dvtemp[0] = tau[0]; */
/*   vdtemp[0] = 1.0; */

/*   /\* read the prior as additional data points *\/ */
/*   if (prior) { */
/*     dcholdc(A0, n_fixed, V); */
/*     for(i = 0; i < n_fixed; i++) { */
/*       X[n_samp+i][n_fixed] = 0; */
/*       for(j = 0; j < n_fixed; j++) { */
/* 	X[n_samp+i][n_fixed] += V[i][j]*beta0[j]; */
/* 	X[n_samp+i][j] = V[i][j]; */
/*       } */
/*     } */
/*   } */

/*   for (j = 0; j < n_random; j++) */
/*     gamma0[j] = 0; */

/*   /\* Gibbs Sampler! *\/ */
/*   for(main_loop = 1; main_loop <= n_gen; main_loop++){ */
/*     /\** STEP 1: Sample Latent Variable **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++){ */
/*       Xbeta[i] = 0; Zgamma[i] = 0; */
/*       for (j = 0; j < n_fixed; j++)  */
/* 	Xbeta[i] += X[i][j]*beta[j];  */
/*       for (j = 0; j < n_random; j++) */
/* 	Zgamma[i] += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j]; */
/*       vitemp[grp[i]]++; */
/*     } */
/*     /\* Sampling tau with MH step *\/ */
/*     if (mh) { */
/*       for (j = 1; j < (n_cat-1); j++)  */
/* 	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1); */
/*       dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000; */
/*       /\* for (j = 0; j < n_cat; j++)  */
/* 	 Rprintf("tau %5d;%14g%14g%14g\n", j, tau[j], dvtemp[j], prop[j]); */
/* 	 Rprintf("\n"); */
/*       *\/ */
/*       for (j = 1; j < (n_cat-1); j++)  */
/* 	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) - */
/* 			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) - */
/* 	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) - */
/* 	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0)); */
/*       for (i = 0; i < n_samp; i++) { */
/* 	if (Y[i] == (n_cat-1))   */
/* 	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-Xbeta[i]-Zgamma[i], 0, 1, 0, 1) - */
/* 	    pnorm(tau[n_cat-2]-Xbeta[i]-Zgamma[i], 0, 1, 0, 1); */
/* 	else if (Y[i] > 0)  */
/* 	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-Xbeta[i]-Zgamma[i], 0, 1, 1, 0) - */
/* 			      pnorm(dvtemp[Y[i]-1]-Xbeta[i]-Zgamma[i], 0, 1, 1, 0)) - */
/* 	    log(pnorm(tau[Y[i]]-Xbeta[i]-Zgamma[i], 0, 1, 1, 0) - */
/* 		pnorm(tau[Y[i]-1]-Xbeta[i]-Zgamma[i], 0, 1, 1, 0)); */
/*       } */
/*       /\* Rprintf("%14g\n", exp(dtemp)); *\/ */
/*       if (unif_rand() < exp(dtemp)) { */
/* 	accept[0]++; */
/* 	for (j = 1; j < n_cat; j++)  */
/* 	  tau[j] = dvtemp[j]; */
/*       } */
/*     }  */

/*     /\* Sampling the Latent Variable *\/ */
/*     if (!mh) { */
/*       Wmin[0] = tau[0]; Wmax[0] = tau[0]-10; */
/*       for (j = 1; j < n_cat; j++) { */
/* 	Wmin[j] = tau[j]; */
/* 	Wmax[j] = tau[j-1]; */
/*       } */
/*     } */

/*     for (i = 0; i < n_samp; i++){ */
/*       if (Y[i] == 0)  */
/* 	W[i] = TruncNorm(Xbeta[i]+Zgamma[i]-1000,0,Xbeta[i]+Zgamma[i],1,0); */
/*       else  */
/* 	W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],Xbeta[i]+Zgamma[i],1,0); */
/*       if (!mh) { */
/* 	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]); */
/* 	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]); */
/*       } */
/*       X[i][n_fixed] = W[i]-Zgamma[i]; */
/*     } */
    
/*     if(!mh) { */
/*       /\* sampling taus without MH-step *\/ */
/*       for (j = 1; j < (n_cat-1); j++)  */
/* 	tau[j] = runif(fmax2(tau[j-1], Wmax[j]),  */
/* 		       fmin2(tau[j+1], Wmin[j+1])); */
/*       tau[n_cat-1] = tau[n_cat-2] + 1000; */
/*     } */

/*     /\** STEP 2: Sample Fixed Effects Given Random Effects **\/ */
/*     bNormalReg(X, beta, vdtemp, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1, */
/* 	       1, 1); */

/*     /\** STEP 3: Update Random Effects Given Fixed Effects **\/ */
/*     for (j = 0; j < n_grp; j++) */
/*       vitemp[j] = 0; */
/*     for (i = 0; i < n_samp; i++) { */
/*       Zgrp[grp[i]][vitemp[grp[i]]][n_random] = W[i]; */
/*       for (j = 0; j < n_fixed; j++)  */
/* 	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j];  */
/*       vitemp[grp[i]]++; */
/*     } */
/*     for (j = 0; j < n_grp; j++) */
/*       bNormalReg(Zgrp[j], gamma[j], vdtemp, vitemp[j], n_random, */
/* 		 1, 1, gamma0, Psi, 0, 0, 1, 1); */

/*     /\** STEP 4: Update Covariance Matrix Given Random Effects **\/ */
/*     for (j = 0; j < n_random; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	mtemp[j][k] = T0[j][k]; */
/*     for (j = 0; j < n_grp; j++) */
/*       for (k = 0; k < n_random; k++) */
/* 	for (l = 0; l < n_random; l++) */
/* 	  mtemp[k][l] += gamma[j][k]*gamma[j][l]; */
/*     dinv(mtemp, n_random, mtemp1); */
/*     rWish(Psi, mtemp1, tau0+n_grp, n_random); */


/*     R_CheckUserInterrupt(); */
/*   } /\* end of Gibbs sampler *\/ */

/*   /\* freeing memory *\/ */
/*   free(gamma0); */
/*   free(Xbeta); */
/*   free(Zgamma); */
/*   FreeMatrix(V, n_fixed); */
/*   free(W); */
/*   free(Wmax); */
/*   free(Wmin); */
/*   FreeMatrix(mtemp, n_random); */
/*   FreeMatrix(mtemp1, n_random); */
/*   free(vdtemp); */
/*   free(dvtemp); */
/*   free(vitemp); */
/* } /\* end of mixed effects ordinal probit *\/ */



/*** 
   A Random Walk Metroplis Sampler for Negative binomial Regression
   with Independent Normal and InvGamma Prior 

   Y_i \sim Negbin(mu_i, sig2) 
      where mu_i = X_i beta is the mean and 
            mu_i + mu_i^2/sig2 is the variance
   prior: 
      beta \sim N(beta_0, A_0^{-1})
      sig2 \sim Gamma(a_0, b_0)
***/

void negbinMetro(int *Y,        /* outcome count variable */
		 double **X,    /* (N x K) covariate matrix */
		 double *beta,  /* K coefficient vector */
		 double *sig2,  /* dispersion parameter */
		 int n_samp,    /* # of obs */
		 int n_cov,     /* # of covariates, K */
		 double *beta0, /* prior mean vector */
		 double **A0,   /* prior precision */
		 double a0,     /* prior shape parameter */
		 double b0,     /* prior scale parameter */
		 double *varb,  /* proposal variances for beta */
		 double vars,   /* proposal variance for sig2 */
		 double *cont,  /* contrast */
		 int n_gen,     /* # of MCMC draws */
		 int *counter,  /* # of acceptance for each parameter
				 */
		 int sig2fixed  /* sig2 fixed? */
		 ) {
  
  int i, j, main_loop;
  double numer, denom;
  double *prop = doubleArray(n_cov);
  double *Xbeta = doubleArray(n_samp);
  double *Xbeta1 = doubleArray(n_samp);

  for (i = 0; i < n_samp; i++) {
    Xbeta[i] = cont[i]; 
    for (j = 0; j < n_cov; j++) 
      Xbeta[i] += X[i][j]*beta[j];
  }
  
  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    /** Sampling beta **/
    for (j = 0; j < n_cov; j++)
      prop[j] = beta[j] + norm_rand()*sqrt(varb[j]);
    /* prior */
    numer = dMVN(prop, beta0, A0, n_cov, 1);
    denom = dMVN(beta, beta0, A0, n_cov, 1);   
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta1[i] = cont[i];
      for (j = 0; j < n_cov; j++) 
	Xbeta1[i] += X[i][j]*prop[j];
      numer += dnegbin(Y[i], exp(Xbeta1[i]), *sig2, 1);
      denom += dnegbin(Y[i], exp(Xbeta[i]), *sig2, 1);
    }
    /* rejection */
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      counter[0]++;
      for (j = 0; j < n_cov; j++)
	beta[j] = prop[j];
      for (i = 0; i < n_samp; i++)
	Xbeta[i] = Xbeta1[i];
    }

    /** Sampling sig2 **/
    if (!sig2fixed) {
      prop[0] = rlnorm(log(sig2[0]), sqrt(vars));
      /* prior */
      numer = dgamma(prop[0], a0, b0, 1);
      denom = dgamma(sig2[0], a0, b0, 1);
      /* likelihood */
      for (i = 0; i < n_samp; i++) {
	numer += dnegbin(Y[i], exp(Xbeta[i]), prop[0], 1);
	denom += dnegbin(Y[i], exp(Xbeta[i]), sig2[0], 1);
      }
      /* proposal distribution */
      denom += dlnorm(prop[0], log(sig2[0]), sqrt(vars), 1);
      numer += dlnorm(sig2[0], log(prop[0]), sqrt(vars), 1);
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	counter[1]++;
	sig2[0] = prop[0];
      }
    }
  }
  
  free(prop);
  free(Xbeta);
  free(Xbeta1);
} /* end of negbinMetro */


void bnegbinMixedMCMC(int *Y,          /* outcome variable */
		      int **Ygrp,      /* outcome variable by group */
		      double **X,      /* model matrix for fixed
					  effects */
		      double ***Zgrp,  /* model matrix for random
					  effects organized by
					  grous */
		      int *grp,        /* group indicator: 0, 1, 2,... */
		      double *beta,    /* fixed effects coefficients */
		      double **gamma,  /* random effects coefficients */
		      double *sig2,    /* dispersion parameter */
		      double **Psi,    /* precision matrix for random
					  effects */
		      int n_samp,      /* # of obs */ 
		      int n_fixed,     /* # of fixed effects */
		      int n_random,    /* # of random effects */
		      int n_grp,       /* # of groups */
		      int max_samp_grp, /* max # of obs per group */
		      double *beta0,   /* prior mean */
		      double **A0,     /* prior precision */
		      double a0,       /* prior shape for sig2 */
		      double b0,       /* prior scale for sig2 */
		      int tau0,        /* prior df for Psi */
		      double **T0,     /* prior scale for Psi */
		      double *varb,    /* proposal variance for beta */
		      double vars,     /* proposal variance for sig2 */
		      double *varg,    /* proposal variance for gamma */
		      int *counter,    /* acceptance counter beta and
					  sig2 2 */
		      int **counterg,  /* acceptance counter for gamma */
		      int n_gen        /* # of gibbs draws */
		      ) {
  
  double *gamma0 = doubleArray(n_random);           /* prior mean for gamma */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  /* storage parameters and loop counters */
  int i, j, k, l, main_loop;  
  int *vitemp = intArray(n_grp);

  /* contrasts */
  double *cont = doubleArray(n_samp);
  double **contM = doubleMatrix(n_grp, max_samp_grp);

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /** STEP 1: Sample Fixed Effects Given Random Effects 
                Also Sample Variance Parameter **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      cont[i] = 0;
      for (j = 0; j < n_random; j++)
	cont[i] += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
      vitemp[grp[i]]++;
    }
    negbinMetro(Y, X, beta, sig2, n_samp, n_fixed, beta0, A0, a0, b0,
		varb, vars, cont, 1, counter, 0);

    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      contM[grp[i]][vitemp[grp[i]]] = 0;
      for (j = 0; j < n_fixed; j++) 
	contM[grp[i]][vitemp[grp[i]]] += X[i][j]*beta[j]; 
      vitemp[grp[i]]++;
    }
    for (j = 0; j < n_grp; j++)
      negbinMetro(Ygrp[j], Zgrp[j], gamma[j], sig2, vitemp[j], n_random,
		  gamma0, Psi, a0, b0, varg, vars, contM[j], 1,
		  counterg[j], 1);

    /** STEP 3: Update Covariance Matrix Given Random Effects **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (j = 0; j < n_grp; j++)
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] += gamma[j][k]*gamma[j][l];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, tau0+n_grp, n_random);

    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(gamma0);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
  free(vitemp);
  free(cont);
  FreeMatrix(contM, n_grp);
} /* end of negative binomial mixed effects model */
