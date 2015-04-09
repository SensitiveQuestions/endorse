#include <string.h>
#include <stdio.h>      
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"
#include "geodist.h"

/* 
   C wrapper code for rtinvchisq() in "endorse" package.
   This function generates a random sample from the truncated inverse
   Chi-squared distribution.
 */
void R2TruncInvChisq(int *n_samp, int *df, double *scale, double *max,
		     double *sample, int *invcdf) {
  int i;

  GetRNGstate();

  for (i = 0; i < *n_samp; i++)
    sample[i] = TruncInvChisq(*df, *scale, *max, *invcdf);

  PutRNGstate();

}



/* 
   C wrapper code for GeoCount() in "endorse" package.
   This function calculate the number of events Y within Z kilometers from points X.
 */
void R2GeoCount (double *dX_lon, double *dX_lat, int *n_X, double *dY_lon, double *dY_lat,
		 int *n_Y, double *Z, int *Store_count) {
  /* loop counters */
  int i, j;
  /* counter for each point */
  int n_temp;
  /* distance */
  double temp;
  /* structs for an event and a point */
  Position X, Y;

  for (i = 0; i < *n_X; i++) {/* for each point of X */
    /* initialize counter */
    n_temp = 0;

    /* location of a point */
    X.lon = dX_lon[i];
    X.lat = dX_lat[i];

    for (j = 0; j < *n_Y; j++) {/* for each event of Y */
      /* location of an event */
      Y.lon = dY_lon[j];
      Y.lat = dY_lat[j];
      
      /* calculate the distance */
      temp = DistanceInMeters(X, Y);
      
      /* update counter */
      if (temp <= *Z * 1000) n_temp++;
    }

    Store_count[i] = n_temp;
  }
}


void R2GeoId (double *dX_lon, double *dX_lat, int *n_X, double *dY_lon, double *dY_lat,
	      int *n_Y, double *Z, int *Store_count) {
  /* loop counters */
  int j;
  /* counter for each point */
  int n_temp = 0;
  /* distance */
  double temp;
  /* structs for an event and a point */
  Position X, Y;

  /* location of a point */
  X.lon = dX_lon[0];
  X.lat = dX_lat[0];
  
  for (j = 0; j < *n_Y; j++) {/* for each event of Y */
    /* location of an event */
    Y.lon = dY_lon[j];
    Y.lat = dY_lat[j];
    
    /* calculate the distance */
    temp = DistanceInMeters(X, Y);
    
    if (temp <= *Z * 1000) {
      Store_count[n_temp] = j+1;
      n_temp++;
    }
  }
}


/* 
   C wrapper code for endorse() in "endorse" pachage.
   This function fits Bayesian measurement model introduced by Bullock, Imai,
   and Shapiro (2011) "Statisitical Analysis of Endorsement Experiment,"
   _Political Analysis_ 19(2): 363-384.
 */
void R2endorse(/*   Data   */
		int *dY,          /* ordinal outcome variable: 0, 1,
				     ..., J-1. Length N * J */
		int *dT,          /* endorsement matrix. Length N * J */
		double *dZ,      /* covariate matirx, length N * M */
		double *dV,     /* village level covariates, length G * P */
		int *village,   /* village indicator for observations, length N. 
				  village[i] \in {0, 1, ..., G}. If not hierarchical,
				  all elements are 0. */
		/*   Data Structure   */
		int *param_indiv,
		int *n_cat,      /* # of categories for each questions:
				    L_j. Length  J*/
		int *param_village,
		/*   Starting values   */
		double *X,      /* vector of ideal points. 
				    No need to include constant
				    Length N */
		double *dS,       /* matrix of support, s_{ij}(k)
				     Length (N * J) */
		double *sbeta,    /* (alpha_j, beta_j), length J * 2 */
		double *stau,     /* cut points for each question;
				     the first cut point is set to 0
				     and the last one set to tau_{L_j - 1}+1000:
				     Length J * max {L_j} */
		double *slambda,  /* lambda, length (J * M) * K  */
		double *somega2,   /* omega2: variance of s given lambda,
				     length J * K */
		double *dtheta,   /* theta, vector of length K * M  */
		double *phi2,     /* phi, vector of length K * M */
		double *kappa,   /* kappa, vector of length P * K */
		double *psi2,    /* psi2, vector of length K */
		double *delta,    /* delta, vector of length M  */
		double *zeta,     /* zeta, vector of length P */
		double *rho2,     /* rho2, vector of length 1 */
		/*  Prior parameters  */
		double *dmu_beta,   /* prior mean of factor loadings */
		double *dA0_beta,     /* prior precision of factor loadings,
				    length 2 * 2: can be set to zero to
				    induce improper prior for beta alone */
		double *dA0_x,     /* known prior precision of ideal points, identical for all i */
		double *dmu_theta, /* prior mean for theta, length K */
		double *dA0_theta, /* prior precision of theta, identical for all k */
		double *dmu_kappa, /* prior mean for kappa, P * K matrix */
		double *dA0_kappa, /* prior precision of kappa, length P * P */
		double *mu_delta,  /* prior mean of delta, length M */
		double *dA0_delta,  /* prior precision of delta, length M * M */
		double *mu_zeta,    /* prior mean of zeta, length P */
		double *dA0_zeta,   /* prior precision of zeta, length P * P */
		double *param_invchisq,
		/* MCMC settings */
		int *settings,
		double *prop,    /* proposal variance for MH step */
		/* Vectors to store sampled parameters */
		double *betaStore,
		double *tauStore,
		double *xStore,
		double *sStore,
		double *lambdaStore,
		double *thetaStore,
		double *kappaStore,
		double *deltaStore,
		double *zetaStore,
		double *omega2Store,
		double *phi2Store,
		double *psi2Store,
		double *sig2Store,
		double *rho2Store,
		/* Vectors to store the last iteration */
		double *betaLast,
		double *tauLast,
		double *xLast,
		double *sLast,
		double *lambdaLast,
		double *thetaLast,
		double *kappaLast,
		double *deltaLast,
		double *zetaLast,
		double *omega2Last,
		double *phi2Last,
		double *psi2Last,
		double *sig2Last,
		double *rho2Last,
		/* Vector to store acceptance ratio */
		double *accept_ratio
		 ){
  /* arguments passed from R */ 
  int n_samp, n_pol, n_dim, n_act, max_n_cat;
  int n_vil, n_vil_dim;
  int n_gen, burn, thin;
  int mda, mh;
  int x_sd, tau_out, s_out, omega2_out, phi2_out, psi2_out, verbose, seed_store, covariates;
  int identical_lambda, hierarchical;

  n_samp = param_indiv[0]; /* # of obs: N */ 
  n_pol = param_indiv[1];     /* # of policies: J */
  n_dim = param_indiv[2]; /* dimension of covariates  */
  n_act = param_indiv[3];    /* # of actors: K */
  max_n_cat = param_indiv[4]; /* max # of categories */
  n_vil = param_village[0];/* number of villages, G. If not hierarchical, 1 */
  n_vil_dim  = param_village[1]; /* dimension of village level covariates, P */
  mda = settings[0]; /* marginal data augmentation? */
  mh = settings[1];         /* Metropolis-Hasting step? */
  x_sd = settings[2];/* Output of x: sd if 1, sample if 0 */
  tau_out = settings[3]; /* Store sampled tau? */
  s_out = settings[4]; /* Store sampled s? */
  omega2_out = settings[5]; /* Store sampled omega2? */
  phi2_out = settings[6]; /* Store sampled phi2? */
  psi2_out = settings[7]; /* Store sampled psi2? */
  verbose = settings[8]; /* Print progress? */
  seed_store = settings[9]; /* Save seed and the last iteration to resume simulation? */
  covariates = settings[10]; /* Covariates for s and x? */
  identical_lambda = settings[11]; /* Identical lambdas across questions? */
  hierarchical = settings[12]; /* hierarchical? */
  n_gen = settings[13];/* # of gibbs draws */
  burn = settings[14];/* # of burnin period */
  thin = settings[15];/* thinning */

  /* loop counters */
  int i, j, k, l, m, n, main_loop, itemp, itempP, progress;
  itempP = (int) ceil( (double) n_gen / 10);
  progress = 1;

  int ibeta = 0, itau = 0, ix = 0, ilambda = 0, itheta = 0, ikappa = 0, idelta = 0, izeta = 0;
  int is = 0, iomega2 = 0, iphi2 = 0, ipsi2 = 0, isig2 = 0, irho2 = 0, keep = 1;
  int ibetaLast = 0, itauLast = 0, ixLast = 0, isLast = 0, ilambdaLast = 0;
  int ithetaLast = 0, ideltaLast = 0, iomega2Last = 0, iphi2Last = 0;
  int ikappaLast = 0;
  double varx;

  /* parameters of the scaled inverse chi-squared prior */
  double s0_omega2, s0_phi2, s0_psi2, s0_sig2, s0_rho2;
  int nu0_omega2, nu0_phi2, nu0_psi2, nu0_sig2, nu0_rho2;
  s0_omega2 = param_invchisq[0];
  nu0_omega2 = param_invchisq[1];
  s0_phi2 = param_invchisq[2];
  nu0_phi2 = param_invchisq[3];
  s0_psi2 = param_invchisq[4];
  nu0_psi2 = param_invchisq[5];
  s0_sig2 = param_invchisq[6];
  nu0_sig2 = param_invchisq[7];
  s0_rho2 = param_invchisq[8];
  nu0_rho2 = param_invchisq[9];

  /* storage */
  double ssr;

  /* storage vectors */
  int *Y_j = intArray(n_samp);
  double *beta = doubleArray(2);
  double *mu_beta = doubleArray(2);
  double *s = doubleArray(1);
  double *var_epsilon = doubleArray(1);
  var_epsilon[0] = 1;
  double *mu_s = doubleArray(1);
  double *x_i = doubleArray(1);
  double *mu_x = doubleArray(1);
  double *lambda_jk = doubleArray(n_dim + n_vil - 1);
  double *omega2_jk = doubleArray(1);
  double *mu_lambda = doubleArray(n_dim + n_vil - 1);
  double *kappa_k = doubleArray(n_vil_dim);
  double *psi2_k = doubleArray(1);
  double *stemp = doubleArray(n_samp * n_pol);
  double *s2temp = doubleArray(n_samp * n_pol);
  double *ztemp = doubleArray(n_samp * n_pol * n_dim);
  int *vtemp = intArray(n_samp * n_pol);
  double *theta_km = doubleArray(n_dim);
  double *phi2_km = doubleArray(1);
  double *mu_theta = doubleArray(1);
  double *sig2_x = doubleArray(1);
  int *accept = intArray(n_pol);
  for (j = 0; j < n_pol; j++)
    accept[j] = 0;
  int *temp_accept = intArray(1);

  /* storage matrices */
  /** data **/
  int **T = intMatrix(n_samp, n_pol); /* treatment */
  int **Y = intMatrix(n_samp, n_pol); /* response */
  double **Z = doubleMatrix(n_samp, n_dim); /* covariates */

  /** starting values **/
  double **S = doubleMatrix(n_samp, n_pol); /* support parameter */
  double **Beta = doubleMatrix(n_pol, 2); /* alpha and beta */
  double **Tau = doubleMatrix(n_pol, max_n_cat); /* cut point */
  double ***Lambda = doubleMatrix3D(n_pol, n_act, n_dim + n_vil - 1); /* lambda */
  double **Omega2 = doubleMatrix(n_pol, n_act); /* omega2 */
  double **Theta = doubleMatrix(n_act, n_dim + n_vil_dim - 1); /* theta */
  double **Kappa = doubleMatrix(n_act * n_pol, n_vil_dim);

  /** prior mean and precision **/
  double **Mu_beta = doubleMatrix(n_pol, 2);
  double **A0_beta = doubleMatrix(2, 2);
  double **A0_s = doubleMatrix(1, 1);
  double **A0_x = doubleMatrix(1, 1);
  double **A0_lambda = doubleMatrix(n_dim + n_vil - 1, n_dim + n_vil - 1);
  double **A0_theta = doubleMatrix(1, 1);
  double **Mu_kappa = doubleMatrix(n_act, n_vil_dim);
  double **A0_kappa = doubleMatrix(n_vil_dim, n_vil_dim);
  double **A0_delta = doubleMatrix(n_dim, n_dim);
  double **A0_zeta = doubleMatrix(n_vil_dim, n_vil_dim);

  /** matrices for regression **/
  double **Ystar = doubleMatrix(n_samp, n_pol); /* latent random utility */
  double **U = doubleMatrix(n_samp+2, 3); /* x_{i} + s_{ij}(k) */
  double **D_s = doubleMatrix(2, 2);
  double **D_x = doubleMatrix(n_pol+1, 2);
  double **D_theta = doubleMatrix(n_pol+1, 2);
  double **D_kappa = doubleMatrix(n_vil + n_vil_dim, n_vil_dim + 1);
  double **D_delta = doubleMatrix(n_samp + n_dim, n_dim + 1);

  int num_D_lambda = 0;
  int num_D_lambda_vil = 0;



  /* get random seed */
  GetRNGstate();

  /* packing data */
  itemp = 0;
  for (j = 0; j < n_pol; j++)
    for (i = 0; i < n_samp; i++)
      Y[i][j] = dY[itemp++];

  itemp = 0;
  for (j = 0; j < n_pol; j++)
    for(i = 0; i < n_samp; i++)
      T[i][j] = dT[itemp++];

  itemp = 0;
  for (m = 0; m < n_dim; m++)
    for (i = 0; i < n_samp; i++)
      Z[i][m] = dZ[itemp++];

  /* packing starting values */
  itemp = 0;
  for (j = 0; j < n_pol; j++)
    for (i = 0; i < n_samp; i++)
      S[i][j] = dS[itemp++];

  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < n_pol; j++)
      Beta[j][m] = sbeta[itemp++];

  itemp = 0;
  for (l = 0; l < max_n_cat; l++)
    for (j = 0; j < n_pol; j++)
      Tau[j][l] = stau[itemp++];

  itemp = 0;
  if (identical_lambda) {

    for (k = 0; k < n_act; k++)
      for (m = 0; m < (n_dim + n_vil - 1); m++)
	Lambda[0][k][m] = slambda[itemp++];

  } else {

    for (k = 0; k < n_act; k++)
      for (j = 0; j < n_pol; j++)
	for (m = 0; m < (n_dim + n_vil - 1); m++)
	  Lambda[j][k][m] = slambda[itemp++];

  }

  itemp = 0;
  if (hierarchical) {
    if (identical_lambda) {
      for (k = 0; k < n_act; k++)
	for (m = 0; m < n_vil_dim; m++)
	  Kappa[k][m] = kappa[itemp++];
    } else {
      for (k = 0; k < n_act * n_pol; k++)
	for (m = 0; m < n_vil_dim; m++)
	  Kappa[k][m] = kappa[itemp++];
    }
  }

  itemp = 0;
  if (identical_lambda) {
    for (k = 0; k < n_act; k++)
      Omega2[0][k] = somega2[itemp++];
  } else {
    for (j = 0; j < n_pol; j++)
      for (k = 0; k < n_act; k++)
	Omega2[j][k] = somega2[itemp++];
  }

  itemp = 0;
  if (hierarchical) {
    for (k = 0; k < n_act; k++)
      for (m = 0; m < (n_dim - 1 + n_vil_dim); m++)
	Theta[k][m] = dtheta[itemp++];
  } else {
    for (k = 0; k < n_act; k++)
      for (m = 0; m < n_dim; m++)
	Theta[k][m] = dtheta[itemp++];
  }

  /* packing prior mean */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < n_pol; j++)
      Mu_beta[j][m] = dmu_beta[itemp++];

  /* packing prior precision */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (n = 0; n < 2; n++)
      A0_beta[n][m] = dA0_beta[itemp++];

  for (m = 0; m < (n_dim + n_vil - 1); m++)
    for (n = 0; n < (n_dim + n_vil - 1); n++)
      A0_lambda[m][n] = 0;

  itemp = 0;
  if (hierarchical) {
    for (m = 0; m < (n_dim - 1); m++)
      for (n = 0; n < (n_dim - 1); n++)
	A0_delta[n][m] = dA0_delta[itemp++];
  } else {
    for (m = 0; m < n_dim; m++)
      for (n = 0; n < n_dim; n++)
	A0_delta[n][m] = dA0_delta[itemp++];
  }

  itemp = 0;
  if (identical_lambda * hierarchical) {
    for (n = 0; n < n_act; n++)
      for (m = 0; m < n_vil_dim; m++)
	Mu_kappa[n][m] = dmu_kappa[itemp++];
  }

  if (hierarchical) {
    itemp = 0;
    for (m = 0; m < n_vil_dim; m++) {
      for (n = 0; n < n_vil_dim; n++) {
	A0_kappa[n][m] = dA0_kappa[itemp];
	A0_zeta[n][m] = dA0_zeta[itemp++];
      }
    }
  }

  A0_x[0][0] = *dA0_x;
  sig2_x[0] = (1 / *dA0_x);


  /* packing covariates */
  if (hierarchical) {
    itemp = 0;
    for (m = 0; m < n_vil_dim; m++)
      for (n = 0; n < n_vil; n++)
	D_kappa[n][m] = dV[itemp++];
  }

  if (covariates) {
    if (hierarchical) {
      for (i = 0; i < n_samp; i++)
	for (m = 0; m < (n_dim - 1); m++)
	  D_delta[i][m] = Z[i][(m + 1)];
    } else {
      for (i = 0; i < n_samp; i++)
	for (m = 0; m < n_dim; m++)
	  D_delta[i][m] = Z[i][m];
    }
  }



  if (identical_lambda) {
    /*** # of observations with treatment k ***/
    for (k = 0; k < n_act; k++) {

      itemp = 0;

      for (j = 0; j < n_pol; j++) {
	for (i = 0; i < n_samp; i++) {
	  if ((k+1) == T[i][j]) {

	    vtemp[itemp] = village[i];
	    itemp++;

	  }
	}
      }

      if (num_D_lambda < itemp) num_D_lambda = itemp;

      if (hierarchical) {
	for (m = 0; m < n_vil; m++) {

	  i = 0;

	  for (l = 0; l < itemp; l++){
	    if (vtemp[l] == m) {

	      i++;

	    }
	  }

	  if (num_D_lambda_vil < i) num_D_lambda_vil = i;

	}
      }
    }
  } else {
    /*** # of observations with treatment k for question j ***/
    for (j = 0; j < n_pol; j++) {
      for (k = 0; k < n_act; k++) {

	itemp = 0;

	for (i = 0; i < n_samp; i++) {
	  if ((k+1) == T[i][j]) {
	    itemp++;
	  }
	}

	if (num_D_lambda < itemp) num_D_lambda = itemp;


	if (hierarchical) {
	  for (m = 0; m < n_vil; m++) {

	    i = 0;

	    for (l = 0; l < itemp; l++){
	      if (vtemp[l] == m) {

		i++;

	      }
	    }

	    if (num_D_lambda_vil < i) num_D_lambda_vil = i;

	  }
	}

      }
    }
  }

  double **D_lambda = doubleMatrix(num_D_lambda + n_dim, n_dim + 1);
  double **D_lambda_vil = doubleMatrix(num_D_lambda_vil + 1, 2);


  /* Gibbs Sampler */
  for (main_loop = 1; main_loop <= n_gen; main_loop++) {
    if (verbose) {
      if (main_loop == 1)
	Rprintf("Start Gibbs Sampler\n");

      if (main_loop  == itempP) {
	Rprintf("%3d percent done.\n   Metropolis acceptance ratios\n",
		progress * 10);
      }
    }

    /** start sampling alpha and beta **/
    for (j = 0; j < n_pol; j++) {

      /*** vectors ***/
      double *tau = doubleArray(n_cat[j]);
      double *MHprop = doubleArray(n_cat[j]-2);

      /*** proposal variance vector ***/
      for (l = 0; l < (n_cat[j]-2); l++)
	MHprop[l] = prop[j];

      /*** response of each question ***/
      for (i = 0; i < n_samp; i++)
	Y_j[i] = Y[i][j];
      
      /*** systematic component of utility ***/
      for (i = 0; i < n_samp; i++) {
	U[i][0] = -1;
	U[i][1] = X[i] + S[i][j];
      }

      /*** starting values of alpha and beta ***/
      for (m = 0; m < 2; m++)
	beta[m] = Beta[j][m];

      /*** starting values of tau ***/
      for (l = 0; l < n_cat[j]; l++)
	tau[l] = Tau[j][l];
      
      /*** prior mean ***/
      for (m = 0; m < 2; m++)
	mu_beta[m] = Mu_beta[j][m];

      /*** set acceptance indicator to 0 ***/
      temp_accept[0] = 0;

      /*** draw posterior ***/
      endorseoprobitMCMC(Y_j, U, beta, tau, n_samp, 2, n_cat[j], 1,
      			  mu_beta, A0_beta, mda, mh, MHprop, temp_accept, 1);

      /*** update acceptance counter ***/
      accept[j] += temp_accept[0];
      accept_ratio[j] = (double) accept[j] / (double) main_loop;

      /*** packing sampled alpha and beta ***/
      for (m = 0; m < 2; m++)
	Beta[j][m] = beta[m];

      /*** packing sampled tau ***/
      for (l = 0; l < n_cat[j]; l++)
	Tau[j][l] = tau[l];

      /*** packing sampled latent random utility ***/
      for (i = 0; i < n_samp; i++)
	Ystar[i][j] = U[i][2];

      /*** storing alpha and beta ***/
      if(main_loop > burn) {
	if(keep == thin) {

	  for (m = 0; m < 2; m++)
	    betaStore[ibeta++] = Beta[j][m];

	  if (tau_out) {
	    for (l = 0; l < (max_n_cat-1); l++)
	      tauStore[itau++] = Tau[j][l];
	  }
	}
      }
      
      /*** store the last iteration ***/
      if (seed_store) {
	if (main_loop == n_gen) {

	  for (m = 0; m < 2; m++)
	    betaLast[ibetaLast++] = Beta[j][m];

	  for (l = 0; l < (max_n_cat - 1); l++)
	    tauLast[itauLast++] = Tau[j][l];

	}
      }


      /** print acceptance ratios  **/
      if (mh * verbose) {
	if (main_loop == itempP)
	  Rprintf("      Cutpoints of question %1d: %4g\n",
		  (j + 1), accept_ratio[j]);
      }

      free(tau);
      free(MHprop);
      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling alpha and beta **/


    /** start sampling s**/
    for (i = 0; i < n_samp; i++){
      for (j = 0; j < n_pol; j++){

    	k = T[i][j];
	
    	/*** if not control, sample s ***/
    	if (k > 0) {

    	  D_s[0][0] = Beta[j][1];
    	  D_s[0][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*X[i];

    	  mu_s[0] = 0;
    	  if (identical_lambda) {
    	    if (hierarchical) {
    	      mu_s[0] += Lambda[0][k - 1][village[i]];

    	      if (covariates) {
    		for (m = 1; m < n_dim; m++)
    		  mu_s[0] += Lambda[0][k - 1][n_vil + m - 1] * Z[i][m];
    	      }

    	    } else {

    	      for (m = 0; m < n_dim; m++)
    		mu_s[0] += Lambda[0][k-1][m] * Z[i][m];

    	    }

    	    A0_s[0][0] = (1 / Omega2[0][k-1]);

    	  } else {
    	    if (hierarchical) {
    	      mu_s[0] += Lambda[j][k - 1][village[i]];

    	      if (covariates) {
    		for (m = 1; m < n_dim; m++)
    		  mu_s[0] += Lambda[j][k - 1][n_vil + m - 1] * Z[i][m];
    	      }

    	    } else {
    	      for (m = 0; m < n_dim; m++)
    		mu_s[0] += Lambda[j][k-1][m] * Z[i][m];
    	    }

    	    A0_s[0][0] = (1 / Omega2[j][k-1]);

    	  }

    	  bNormalReg(D_s, s, var_epsilon, 1, 1, 1, 1, mu_s, A0_s,
    	  	     1, 1, 1, 1, 0);

    	  /*** packing sampled s ***/
    	  S[i][j] = s[0];
    	}

    	/*** storing sampled s  ***/
    	if (s_out) {
    	  if(main_loop > burn) {
    	    if(keep == thin) {
    	      sStore[is++] = S[i][j];
    	    }
    	  }
    	}
    	/*** storing the last iteration ***/
    	if (seed_store) {
    	  if (main_loop == n_gen) {
    	    sLast[isLast++] = S[i][j];
    	  }
    	}

    	R_FlushConsole();
    	R_CheckUserInterrupt();
      }
    }/** end of sampling s **/


    /** start sampling x **/
    for (i = 0; i < n_samp; i++) {

      for (j = 0; j < n_pol; j++) {
    	D_x[j][0] = Beta[j][1];
    	D_x[j][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1] * S[i][j];
      }
      
      /*** prior mean ***/
      if (hierarchical) {

    	mu_x[0] = delta[village[i]];

    	if (covariates) {
    	  for (m = 0; m < (n_dim - 1); m++)
    	    mu_x[0] += Z[i][(m + 1)] * delta[n_vil + m];
    	}

      }	else if (covariates) {
    	mu_x[0] = 0;
    	for (m = 0; m < n_dim; m++)
    	  mu_x[0] += Z[i][m] * delta[m];
      } else {
    	mu_x[0] = mu_delta[i];
      }
      
      bNormalReg(D_x, x_i, var_epsilon, n_pol, 1, 1, 1, mu_x, A0_x,
      		 1, 1, 1, 1, 0);

      /*** packing sampled x ***/
      X[i] = x_i[0];

      R_FlushConsole();
      R_CheckUserInterrupt();
    }

    /*** storing sampled x ***/
    if(main_loop > burn) {
      if (keep == thin) {
    	if (x_sd){
    	  varx = var(X, n_samp, 1);
    	  xStore[ix++] = sqrt(varx);
    	} else {
    	  for (i = 0; i < n_samp; i++)
    	    xStore[ix++] = X[i];
    	}
      }
    }
    /**** if the output may be used for updating ****/
    if (seed_store) {
    if (main_loop == n_gen) {
    	for (i = 0; i < n_samp; i++)
    	  xLast[ixLast++] = X[i];
      }
    } /** end of sampling x **/



    if (identical_lambda) { /** start sampling lambda,
    				 lambdas are identical across policies **/
      for (k = 0; k < n_act; k++) {
    	/*** # of observations with treatment k for question j ***/
    	n = 0;
    	itemp = 0;
    	for (j = 0; j < n_pol; j++) {
    	  for (i = 0; i < n_samp; i++) {
    	    if ((k+1) == T[i][j]) {
    	      stemp[n] = S[i][j];

    	      /* individual level covariates except for the intercept */
    	      for (m = 1; m < n_dim; m++)
    		ztemp[itemp++] = Z[i][m];

    	      /* village specific intercept */
    	      vtemp[n] = village[i];

    	      n++;
    	    }
    	  }
    	}
	
    	if (n == 0)
    	  continue;


    	if (hierarchical) {
	  
    	  if (covariates) {
    	    /* sampling coefficients */
    	    itemp = 0;
    	    for (l = 0; l < n; l++) {

    	      /* packing covariates except for the intercept */
    	      for (m = 0; m < (n_dim - 1); m++)
    		D_lambda[l][m] = ztemp[itemp++];

    	      /* regressand, s */
    	      D_lambda[l][n_dim - 1] = stemp[l] - Lambda[0][k][vtemp[l]];
    	    }

    	    /* prior */
    	    for (m = 0; m < (n_dim - 1); m++) {
    	      mu_lambda[m] = dmu_theta[(k * (n_dim - 1) + m)];
    	      A0_lambda[m][m] = dA0_theta[m];
    	    }

    	    omega2_jk[0] = Omega2[0][k];

    	    bNormalReg(D_lambda, lambda_jk, omega2_jk, n, (n_dim - 1), 1, 1,
    		       mu_lambda, A0_lambda, 1, s0_omega2, nu0_omega2, 0, 0);

    	    Omega2[0][k] = omega2_jk[0];

    	    for (m = 0; m < (n_dim - 1); m++)
    	      Lambda[0][k][n_vil + m] = lambda_jk[m];
    	  } else {
    	    /* sampling omega2 when hierarchical model without individual level covariate */
    	    ssr = 0;
    	    for (l = 0; l < n; l++) {
    	      s2temp[l] = (stemp[l] - Lambda[0][k][vtemp[l]]) * (stemp[l] - Lambda[0][k][vtemp[l]]);
    	      ssr += s2temp[l];
    	    }
    	    Omega2[0][k] = (ssr + nu0_omega2 * s0_omega2) / rchisq((double)n+nu0_omega2);
    	  }

    	  itemp = 0;
    	  for (l = 0; l < n; l++)
    	    for (m = 0; m < (n_dim - 1); m++)
    	      D_lambda[l][m] = ztemp[itemp++];
	  

    	  /* sampling village specific intercept */
    	  for (m = 0; m < n_vil; m++) {
    	    itemp = 0;
    	    for (l = 0; l < n; l++){
    	      if (vtemp[l] == m) {

    		D_lambda_vil[itemp][1] = stemp[l];

    		if (covariates) {
    		  for (j = 0; j < (n_dim - 1); j++) {
    		    D_lambda_vil[itemp][1] -= D_lambda[l][j] * Lambda[0][k][n_vil + j];
    		  }
    		}

    		D_lambda_vil[itemp][0] = 1;

    		itemp++;
    	      }
    	    }


    	    /** add prior **/
    	    mu_lambda[0] = 0;
    	    for (j = 0; j < n_vil_dim; j++)
    	      mu_lambda[0] += D_kappa[m][j] * Kappa[k][j];

    	    A0_lambda[0][0] = (1 / psi2[k]);

    	    omega2_jk[0] = Omega2[0][k];

    	    bNormalReg(D_lambda_vil, lambda_jk, omega2_jk, itemp, 1, 1, 1,
    		       mu_lambda, A0_lambda, 0, 1, 1, 1, 0);

    	    Lambda[0][k][m] = lambda_jk[0];
    	  }

    	} else { /* non-hierarchical model */
    	  itemp = 0;
    	  for (l = 0; l < n; l++) {
    	    /* packing covariates except for the intercept */
    	    for (m = 0; m < (n_dim - 1); m++)
    	      D_lambda[l][n_vil + m] = ztemp[itemp++];

    	    /* intercept */
    	    D_lambda[l][0] = 1;

    	    /* regressand, s */
    	    D_lambda[l][n_dim + n_vil - 1] = stemp[l];
    	  }
	  
    	  /* prior */
    	  for (m = 0; m < n_dim; m++) {
    	    mu_lambda[m] = dmu_theta[m];
    	    A0_lambda[m][m] = dA0_theta[m];
    	  }
	
    	  omega2_jk[0] = Omega2[0][k];

    	  bNormalReg(D_lambda, lambda_jk, omega2_jk, n, (n_dim + n_vil - 1), 1, 1,
    		     mu_lambda, A0_lambda, 0, s0_omega2, nu0_omega2, 0, 0);

    	  for (m = 0; m < (n_dim + n_vil - 1); m++)
    	    Lambda[0][k][m] = lambda_jk[m];

    	  Omega2[0][k] = omega2_jk[0];
    	}


    	/*** storing sampled lambda ***/
    	if (main_loop > burn) {
    	  if (keep == thin) {
    	    if (hierarchical) {
    	      for (m = 0; m < (n_dim + n_vil - 1); m++)
    		lambdaStore[ilambda++] = Lambda[0][k][m];
    	    } else {
    	      for (m = 0; m < n_dim; m++)
    		lambdaStore[ilambda++] = Lambda[0][k][m];
    	    }

    	    if (omega2_out) omega2Store[iomega2++] = Omega2[0][k];
    	  }
    	}

    	if (seed_store) {
    	  if (main_loop == n_gen) {
    	    for (m = 0; m < (n_dim + n_vil - 1); m++)
    	      lambdaLast[ilambdaLast++] = Lambda[0][k][m];

    	    omega2Last[iomega2Last++] = Omega2[0][k];
    	  }
    	}

    	R_FlushConsole();
    	R_CheckUserInterrupt();
      }
    } else { /** start sampling lambda and omega2,
    		 lambdas are NOT identical across policies **/
      for (k = 0; k < n_act; k++) {
    	for (j = 0; j < n_pol; j++) {

    	  /*** # of observations with treatment k for question j ***/
    	  n = 0;
    	  itemp = 0;
    	  for (i = 0; i < n_samp; i++) {
    	    if ((k+1) == T[i][j]) {
    	      stemp[n] = S[i][j];

    	      for (m = 0; m < n_dim; m++)
    		ztemp[itemp++] = Z[i][m];

    	      /* village specific intercept */
    	      vtemp[n] = village[i];

    	      n++;
    	    }
    	  }
	
    	  if (n == 0)
    	    continue;

    	  if (hierarchical) {
	  
    	    if (covariates) {
    	      /* sampling coefficients */
    	      itemp = 0;
    	      for (l = 0; l < n; l++) {

    		/* packing covariates except for the intercept */
    		for (m = 0; m < (n_dim - 1); m++)
    		  D_lambda[l][m] = ztemp[itemp++];

    		/* regressand, s */
    		D_lambda[l][n_dim - 1] = stemp[l] - Lambda[j][k][vtemp[l]];
    	      }

    	      /* prior */
    	      for (m = 0; m < (n_dim - 1); m++) {
    		mu_lambda[m] = Theta[k][(n_vil_dim + n_dim - 1) + m];
    		A0_lambda[m][m] = 1/phi2[k * (n_vil_dim + n_dim - 1) + m];
    	      }

    	      omega2_jk[0] = Omega2[j][k];

    	      bNormalReg(D_lambda, lambda_jk, omega2_jk, n, (n_dim - 1), 1, 1,
    			 mu_lambda, A0_lambda, 1, s0_omega2, nu0_omega2, 0, 0);

    	      Omega2[j][k] = omega2_jk[0];

    	      for (m = 0; m < (n_dim - 1); m++)
    		Lambda[j][k][n_vil + m] = lambda_jk[m];
    	    } else {
    	      /* sampling omega2 when hierarchical model without individual level covariate */
    	      ssr = 0;
    	      for (l = 0; l < n; l++) {
    		s2temp[l] = (stemp[l] - Lambda[j][k][vtemp[l]]) * (stemp[l] - Lambda[j][k][vtemp[l]]);
    		ssr += s2temp[l];
    	      }
    	      Omega2[j][k] = (ssr + nu0_omega2 * s0_omega2) / rchisq((double)n+nu0_omega2);
    	    }

    	    itemp = 0;
    	    for (l = 0; l < n; l++)
    	      for (m = 0; m < (n_dim - 1); m++)
    		D_lambda[l][m] = ztemp[itemp++];



    	    /* sampling village specific intercept */
    	    for (m = 0; m < n_vil; m++) {
    	      itemp = 0;
    	      for (l = 0; l < n; l++){
    		if (vtemp[l] == m) {

    		  D_lambda_vil[itemp][1] = stemp[l];

    		  if (covariates) {
    		    for (i = 0; i < (n_dim - 1); i++) {
    		      D_lambda_vil[itemp][1] -= D_lambda[l][i] * Lambda[j][k][n_vil + i];
    		    }
    		  }

    		  D_lambda_vil[itemp][0] = 1;

    		  itemp++;
    		}
    	      }



    	      /** add prior **/
    	      mu_lambda[0] = 0;
    	      for (i = 0; i < n_vil_dim; i++)
    		mu_lambda[0] += D_kappa[m][i] * Kappa[k * n_pol + j][i];

    	      A0_lambda[0][0] = (1 / psi2[k * n_pol + j]);

    	      omega2_jk[0] = Omega2[j][k];

    	      bNormalReg(D_lambda_vil, lambda_jk, omega2_jk, itemp, 1, 1, 1,
    	      		 mu_lambda, A0_lambda, 0, 1, 1, 1, 0);

    	      Lambda[j][k][m] = lambda_jk[0];

    	    }


    	  } else {
	
    	    itemp = 0;
    	    for (l = 0; l < n; l++) {

    	      for (m = 0; m < n_dim; m++)
    		D_lambda[l][m] = ztemp[itemp++];
	      
    	      D_lambda[l][n_dim] = stemp[l];
    	    }

    	    for (m = 0; m < n_dim; m++) {
    	      mu_lambda[m] = Theta[k][m];
    	      A0_lambda[m][m] = (1 / phi2[(k * n_dim) + m]);
    	    }

    	    for (m = 0; m < n_dim; m++)
    	      lambda_jk[m] = Lambda[j][k][m];

    	    omega2_jk[0] = Omega2[j][k];
	
    	    bNormalReg(D_lambda, lambda_jk, omega2_jk, n, n_dim, 1, 1, mu_lambda,
    		       A0_lambda, 0, s0_omega2, nu0_omega2, 0, 0);

    	    for (m = 0; m < n_dim; m++)
    	      Lambda[j][k][m] = lambda_jk[m];

    	    Omega2[j][k] = omega2_jk[0];

    	  }

    	  if (main_loop > burn) {
    	    if (keep == thin) {
    	      for (m = 0; m < (n_dim + n_vil - 1); m++)
    		lambdaStore[ilambda++] = Lambda[j][k][m];

    	      if (omega2_out) omega2Store[iomega2++] = Omega2[j][k];
    	    }
    	  }

      /**** if the chain may be updated later ****/
    	  if (seed_store) {
    	    if (main_loop == n_gen) {
    	      for (m = 0; m < n_dim + n_vil - 1; m++)
    		lambdaLast[ilambdaLast++] = Lambda[j][k][m];

    	      omega2Last[iomega2Last++] = Omega2[j][k];
    	    }
    	  }

    	  R_FlushConsole();
    	  R_CheckUserInterrupt();
    	}
      }
    }/** end of sampling lambda and omega2 **/




    if (!identical_lambda) {/** start sampling theta and phi2 **/
      for (k = 0; k < n_act; k++) {

    	if (hierarchical) {
    	  /** individual level parameters **/
    	  if (covariates) {
    	    for (m = 0; m < n_dim - 1; m++) {

    	      for (j = 0; j < n_pol; j++) {
    		D_theta[j][0] = 1;
    		D_theta[j][1] = Lambda[j][k][n_vil + m];
    	      }
	
    	      /*** prior mean for theta_{k} ***/
    	      mu_theta[0] = dmu_theta[k * (n_dim - 1 + n_vil_dim) + m];

    	      A0_theta[0][0] = dA0_theta[m];

    	      theta_km[0] = Theta[k][m];
    	      phi2_km[0] = phi2[k * (n_dim - 1 + n_vil_dim) + m];

    	      bNormalReg(D_theta, theta_km, phi2_km, n_pol, 1, 1, 1,
    			 mu_theta, A0_theta, 0, s0_phi2, nu0_phi2, 0, 0);

    	      /*** packing sampled theta ***/
    	      Theta[k][m] = theta_km[0];

    	      /*** packing sampled phi2 ***/
    	      phi2[k * (n_dim - 1 + n_vil_dim) + m] = phi2_km[0];

    	      /*** storing sampled theta ***/
    	      if (main_loop > burn) {
    		if (keep == thin) {
    		  thetaStore[itheta++] = Theta[k][m];

    		  if (phi2_out) phi2Store[iphi2++] = phi2[k * (n_dim - 1 + n_vil_dim) + m];
    		}
    	      }

    	      if (seed_store) {
    		if (main_loop == n_gen) {
    		  thetaLast[ithetaLast++] = Theta[k][m];

    		  phi2Last[iphi2Last++] = phi2[k * (n_dim - 1 + n_vil_dim) + m];
    		}
    	      }
    	    }
    	  }

    	  /** group level parameters **/
    	  for (m = 0; m < n_vil_dim; m++) {

    	    for (j = 0; j < n_pol; j++) {
    	      D_theta[j][0] = 1;
    	      D_theta[j][1] = Kappa[k * n_pol + j][m];
    	    }
	
    	    /*** prior mean for theta_{k} ***/
    	    mu_theta[0] = dmu_theta[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m];

    	    A0_theta[0][0] = dA0_theta[n_dim - 1 + m];

    	    theta_km[0] = Theta[k][n_dim - 1 + m];
    	    phi2_km[0] = phi2[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m];

    	    bNormalReg(D_theta, theta_km, phi2_km, n_pol, 1, 1, 1,
    		       mu_theta, A0_theta, 0, s0_phi2, nu0_phi2, 0, 0);

    	    /*** packing sampled theta ***/
    	    Theta[k][n_dim - 1 + m] = theta_km[0];

    	    /*** packing sampled phi2 ***/
    	    phi2[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m] = phi2_km[0];

    	    /*** storing sampled theta ***/
    	    if (main_loop > burn) {
    	      if (keep == thin) {
    		thetaStore[itheta++] = Theta[k][n_dim - 1 + m];

    		if (phi2_out) phi2Store[iphi2++] = phi2[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m];
    	      }
    	    }

    	    if (seed_store) {
    	      if (main_loop == n_gen) {
    		thetaLast[ithetaLast++] = Theta[k][n_dim - 1 + m];

    		phi2Last[iphi2Last++] = phi2[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m];
    	      }
    	    }
    	  }

    	  R_FlushConsole();
    	  R_CheckUserInterrupt();
	  
    	} else {
    	  for (m = 0; m < n_dim; m++) {

    	    for (j = 0; j < n_pol; j++) {
    	      D_theta[j][0] = 1;
    	      D_theta[j][1] = Lambda[j][k][m];
    	    }
	
    	    /*** prior mean for theta_{k} ***/
    	    mu_theta[0] = dmu_theta[m];

    	    A0_theta[0][0] = dA0_theta[m];

    	    theta_km[0] = Theta[k][m];
    	    phi2_km[0] = phi2[(k * n_dim) + m];

    	    bNormalReg(D_theta, theta_km, phi2_km, n_pol, 1, 1, 1,
    		       mu_theta, A0_theta, 0, s0_phi2, nu0_phi2, 0, 0);

    	    /*** packing sampled theta ***/
    	    Theta[k][m] = theta_km[0];

    	    /*** packing sampled phi2 ***/
    	    phi2[(k * n_dim) + m] = phi2_km[0];

    	    /*** storing sampled theta ***/
    	    if (main_loop > burn) {
    	      if (keep == thin) {
    		thetaStore[itheta++] = Theta[k][m];

    		if (phi2_out) phi2Store[iphi2++] = phi2[(k * n_dim) + m];
    	      }
    	    }

    	    if (seed_store) {
    	      if (main_loop == n_gen) {
    		thetaLast[ithetaLast++] = Theta[k][m];

    		phi2Last[iphi2Last++] = phi2[(k * n_dim) + m];
    	      }
    	    }
    	  }

    	  R_FlushConsole();
    	  R_CheckUserInterrupt();
    	}
      }
    }/** end of sampling theta and phi2 **/


    if (hierarchical) {/** start sampling kappa and psi2, village level regression */
      for (k = 0; k < n_act; k++) {

    	if (identical_lambda) {

    	  l = 1;

    	} else {

    	  l = n_pol;
	  
	  for (m = 0; m < n_vil_dim; m++)
	    for (n = 0; n < n_vil_dim; n++)
	      A0_kappa[n][m] = 0;

	  for (m = 0; m < n_vil_dim; m++) {
	    Mu_kappa[k][m] = Theta[k][n_dim - 1 + m];
	    A0_kappa[m][m] = phi2[k * (n_dim - 1 + n_vil_dim) + n_dim - 1 + m];
	  }

    	}

    	for (j = 0; j < l; j++) {

    	  for (n = 0; n < n_vil; n++)
    	    D_kappa[n][n_vil_dim] = Lambda[j][k][n];

    	  psi2_k[0] = psi2[k * l + j];

    	  bNormalReg(D_kappa, kappa_k, psi2_k, n_vil, n_vil_dim, 1, 1, Mu_kappa[k],
    	  	     A0_kappa, 1, s0_psi2, nu0_psi2, 0, 0);

    	  /* re-packing covariates */
    	  itemp = 0;
    	  for (m = 0; m < n_vil_dim; m++)
    	    for (n = 0; n < n_vil; n++)
    	      D_kappa[n][m] = dV[itemp++];

    	  for (m = 0; m < n_vil_dim; m++)
    	    Kappa[k * l + j][m] = kappa_k[m];

    	  psi2[k * l + j] = psi2_k[0];

    	  if (main_loop > burn) {
    	    if (keep == thin) {

    	      for (m = 0; m < n_vil_dim; m++)
    	      	kappaStore[ikappa++] = Kappa[k * l + j][m];

    	      if (psi2_out) psi2Store[ipsi2++] = psi2[k * l + j];
    	    }
    	  }

    	  if (seed_store) {
    	    if (main_loop == n_gen) {
    	      for (m = 0; m < n_vil_dim; m++)
    		kappaLast[ikappaLast++] = Kappa[k * l + j][m];

    	      psi2Last[k] = psi2[k * l + j];
    	    }
    	  }
    	}
      }
    }


    /** Start sampling delta **/
    if (hierarchical) { /** hierarchical model **/

      if (covariates) {
    	/**  packing covariates **/
    	for (i = 0; i < n_samp; i++)
    	  for (m = 0; m < (n_dim - 1); m++)
    	    D_delta[i][m] = Z[i][(m + 1)];

    	/** coefficients on individual level covariates **/
    	for (i = 0; i < n_samp; i++)
    	  D_delta[i][(n_dim - 1)] = X[i] - delta[village[i]];

    	bNormalReg(D_delta, lambda_jk, sig2_x, n_samp, (n_dim - 1), 1, 1,
    		   mu_delta, A0_delta, 1, s0_sig2, nu0_sig2, 0, 0);

    	for (m = 0; m < (n_dim - 1); m++)
    	  delta[n_vil + m] = lambda_jk[m];
      } else {
    	/** sampling sig2 if hierarchical model without covariates **/
    	ssr = 0;
    	for (i = 0; i < n_samp; i++) {
    	  s2temp[i] = (X[i] - delta[village[i]]) * (X[i] - delta[village[i]]);
    	  ssr += s2temp[i];
    	}
    	sig2_x[0] = (ssr + s0_sig2 * nu0_sig2) / rchisq((double)n_samp+nu0_sig2);
      }

      /** random intercept for each village **/
      for (m = 0; m < n_vil; m++) {
    	itemp = 0;
    	for (i = 0; i < n_samp; i++) {
    	  if (village[i] == m) {
    	    D_delta[itemp][1] = X[i];

    	    if (covariates) {
    	      for (l = 0; l < (n_dim - 1); l++)
    		D_delta[itemp][1] -= Z[i][(l + 1)] * delta[n_vil + l];
    	    }

    	    D_delta[itemp][0] = 1;
	      
    	    itemp++;
    	  }
    	}


    	/** prior for village specific intercept **/
    	mu_theta[0] = 0;
    	for (l = 0; l < n_vil_dim; l++)
    	  mu_theta[0] += D_kappa[m][l] * zeta[l];

    	A0_theta[0][0] = (1 / rho2[0]);
	  
    	/** regression **/
    	bNormalReg(D_delta, lambda_jk, sig2_x, itemp, 1, 1, 1,
    		   mu_theta, A0_theta, 1, 1, 1, 1, 0);

    	delta[m] = lambda_jk[0];
	  
      }

      if (main_loop > burn) {
    	if (keep == thin) {
    	  for (m = 0; m < (n_vil + n_dim - 1); m++)
    	    deltaStore[idelta++] = delta[m];

    	  sig2Store[isig2++] = sig2_x[0];
    	}
      }

      if (seed_store) {
    	if (main_loop == n_gen) {
    	  for (m = 0; m < (n_vil + n_dim - 1); m++)
    	    deltaLast[ideltaLast++] = delta[m];

    	  sig2Last[0] = sig2_x[0];
    	}
      }


    } else if (covariates) { /** non-hierarchical model **/

      for (i = 0; i < n_samp; i++)
    	D_delta[i][n_dim] = X[i];

      bNormalReg(D_delta, delta, sig2_x, n_samp, n_dim, 1, 1,
    		 mu_delta, A0_delta, 1, 1, 1, 1, 0);

      if (main_loop > burn) {
    	if (keep == thin) {
    	  for (m = 0; m < n_dim; m++)
    	    deltaStore[idelta++] = delta[m];
    	}
      }

      if (seed_store) {
    	if (main_loop == n_gen) {
    	  for (m = 0; m < n_dim; m++)
    	    deltaLast[ideltaLast++] = delta[m];
    	}
      }

    }
    /** end of sampling delta **/

    if (hierarchical) { /** sampling zeta **/

      for (m = 0; m < n_vil; m++)
      	D_kappa[m][n_vil_dim] = delta[m];

      bNormalReg(D_kappa, zeta, rho2, n_vil, n_vil_dim, 1, 1, mu_zeta, A0_zeta,
      		 1, s0_rho2, nu0_rho2, 0, 0);

      	/* re-packing covariates */
    	itemp = 0;
    	for (m = 0; m < n_vil_dim; m++)
    	  for (n = 0; n < n_vil; n++)
    	    D_kappa[n][m] = dV[itemp++];
	
    	if (main_loop > burn) {
    	  if (keep == thin) {
    	    for (m = 0; m < n_vil_dim; m++)
    	      zetaStore[izeta++] = zeta[m];

    	    rho2Store[irho2++] = rho2[0];
    	  }
    	}

    	if (seed_store) {
    	  if (main_loop == n_gen) {
    	    for (m = 0; m < n_vil_dim; m++)
    	      zetaLast[m] = zeta[m];

    	    rho2Last[0] = rho2[0];
    	  }
    	}
    }


    /** update thinning counter **/
    if (keep == thin) {
      keep = 1;
    } else {
      keep++;
    }

    /** update printing counter **/
    if (verbose) {
      if (main_loop == itempP) {
	progress += 1;
	itempP = (int) ceil( (double) progress * n_gen / 10);      
      }
    }

    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  if (verbose)
    Rprintf("End Gibbs Sampler\n");

  PutRNGstate();

  /* freeing memory */
  free(Y_j);
  free(beta);
  free(mu_beta);
  free(s);
  free(var_epsilon);
  free(mu_s);
  free(x_i);
  free(mu_x);
  free(lambda_jk);
  free(omega2_jk);
  free(mu_lambda);
  free(kappa_k);
  free(psi2_k);
  free(stemp);
  free(ztemp);
  free(vtemp);
  free(theta_km);
  free(phi2_km);
  free(mu_theta);
  free(sig2_x);
  free(accept);
  free(temp_accept);

  FreeintMatrix(T, n_samp);
  FreeintMatrix(Y, n_samp);
  FreeMatrix(Z, n_samp);
  FreeMatrix(S, n_samp);
  FreeMatrix(Beta, n_pol);
  FreeMatrix(Tau, n_pol);
  FreeMatrix(Theta, n_act);
  Free3DMatrix(Lambda, n_pol, n_act);
  FreeMatrix(Omega2, n_pol);
  FreeMatrix(Kappa, n_act);
  FreeMatrix(Mu_beta, n_pol);
  FreeMatrix(A0_beta, 2);
  FreeMatrix(Ystar, n_samp);
  FreeMatrix(U, n_samp+2);
  FreeMatrix(D_s, 2);
  FreeMatrix(A0_s, 1);
  FreeMatrix(D_x, n_pol+1);
  FreeMatrix(A0_x, 1);
  FreeMatrix(A0_lambda, n_dim + n_vil - 1);
  FreeMatrix(D_theta, n_pol+1);
  FreeMatrix(A0_theta, 1);
  FreeMatrix(A0_delta, n_dim);
  FreeMatrix(D_delta, n_samp + n_dim);
  FreeMatrix(D_kappa, n_vil + n_vil_dim);
  FreeMatrix(D_lambda, num_D_lambda + n_dim);
  FreeMatrix(D_lambda_vil, num_D_lambda_vil + 1);
  FreeMatrix(Mu_kappa, n_act);
  FreeMatrix(A0_kappa, n_vil_dim);
  FreeMatrix(A0_zeta, n_vil_dim);
}
