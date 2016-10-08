\name{endorse}
\alias{endorse}
\title{
  Fitting the Measurement Model of Political Support via Markov Chain
  Monte Carlo
}
\description{
  This function generates a sample from the posterior distribution of the
  measurement model of political support. Individual-level covariates
  may be included in the model. The details of the model are given under
  `Details'. See also Bullock et al. (2011).
}
\usage{
endorse(Y, data, data.village = NA, village = NA, treat = NA,
        na.strings = 99, identical.lambda = TRUE, covariates = FALSE,
        formula.indiv = NA, hierarchical = FALSE, formula.village = NA,
        x.start = 0, s.start = 0, beta.start = 1, tau.start = NA,
        lambda.start = 0, omega2.start = .1, theta.start = 0,
        phi2.start = .1, kappa.start = 0, psi2.start = 1,
        delta.start = 0, zeta.start = 0, rho2.start = 1, mu.beta = 0,
        mu.x = 0, mu.theta = 0, mu.kappa = 0, mu.delta = 0, mu.zeta = 0,
        precision.beta = 0.04, precision.x = 1, precision.theta = 0.04,
        precision.kappa = 0.04, precision.delta = 0.04,
        precision.zeta = 0.04, s0.omega2= 1, nu0.omega2 = 10,
        s0.phi2 = 1, nu0.phi2 = 10, s0.psi2 = 1, nu0.psi2 = 10,
        s0.sig2 = 1, nu0.sig2 = 400, s0.rho2 = 1, nu0.rho2 = 10,
        MCMC = 20000, burn = 1000, thin = 1,  mh = TRUE, prop = 0.001,
        x.sd = TRUE, tau.out = FALSE, s.out = FALSE, omega2.out = TRUE,
        phi2.out = TRUE, psi2.out = TRUE, verbose = TRUE,
        seed.store = FALSE, update = FALSE, update.start = NULL)
}

\arguments{

  \item{Y}{a list of the variable names for the responses. It should
    take the following form:
    \code{list(Q1 = c("varnameQ1.1", ...), Q2 = c("varnameQ2.1", ...), ...)}.
    If \code{treat} is \code{NA}, the first
    variable for each question should be the responses of the
    control observations while each of the other variables should
    correspond to each endorser. \code{treat} should be supplied if
    only one variable name is provided for a question in this argument.
  }

  \item{data}{data frame containing the individual-level variables.
    The cases must be complete, i.e., no \code{NA}'s are allowed.
  }

  \item{data.village}{data frame containing the village-level variables.
    The cases must be complete, i.e., no \code{NA}'s are allowed. If 
    auxiliary information is included, the data frame should include only
    a unique identifier for each unit for which prediction is desired.
  }

  \item{village}{character.  The variable name of the village indicator
    in the individual-level data. If auxiliary information is included,
    this should correspond to the variable name of the units at which
    prediction is desired. 
  }

  \item{treat}{An optional matrix of non negative integers indicating
    the treatment
    status of each observation and each question.
    Rows are observations
    and columns are questions. 0 represents the control status while
    positive integers indicate treatment statuses.
    If \code{treat} is set to \code{NA}, the function generates the    
    treatment matrix using \code{Y}.
    The default is \code{NA}.
  }
  
  \item{na.strings}{a scalar or a vector indicating the values of the
    response variable that are to be interpreted as ``Don't Know'' or
    ``Refused to Answer.'' The value should not be \code{NA} unless
    \code{treat} is provided, because \code{NA}'s are interpreted as the
    response to the question with another endorsement. Default is
    \code{99}.
  }

  \item{identical.lambda}{logical.  If \code{TRUE}, the model with a common
    lambda across questions will be fitted.  The default is \code{TRUE}.
  }

  \item{covariates}{logical. If \code{TRUE}, the model includes
    individual-level covariates. The default is \code{FALSE}.
  }

  \item{formula.indiv}{a symbolic description specifying the individual level 
    covariates for the support parameter and the ideal points. The formula
    should be one-sided, e.g. \code{ ~ Z1 + Z2}.
  }

  \item{hierarchical}{logical.  IF \code{TRUE}, the hierarchical model with
    village level predictors will be fitted.  The default is \code{FALSE}.
  }

  \item{formula.village}{a symbolic description specifying the village level 
    covariates for the support parameter and the ideal points. The formula
    should be one-sided, e.g. \code{ ~ V1 + V2}.
  }
  
  \item{h}{Auxiliary data functionality. Optional named numeric vector with 
  length equal to number of groups. Names correspond to group labels and 
  values correspond to auxiliary moments.
  }
  
  \item{group}{Auxiliary data functionality. Optional character vector of 
  group labels with length equal to number of observations.
  }

  \item{x.start}{starting values for the ideal points vector \eqn{x}. If
    \code{x.start} is set to a scalar, the starting values for the ideal
    points of all respondents will be set to the scalar. If
    \code{x.start} is a vector of the same length as the number of
    observations, then this vector will be used as the starting
    values. The default is \code{0}.
  }

  \item{s.start}{starting values for the support parameter, \eqn{s_ijk}.
    If \code{s.start} is set to a scalar, the starting values for the
    support parameter of all respondents and all questions will be the
    scalar.
    If \code{s.start} is set to a matrix, it should have the same number
    of rows as the number of observations and the same number of columns
    as the number of questions. Also, the value should be zero for the
    control condition.
    The default is \code{0}.
  }

  \item{beta.start}{starting values for the question related parameters,
    \eqn{\alpha_j} and \eqn{\beta_j}.
    If \code{beta.start} is set to a scalar, the starting values for the
    support parameter of all respondents and all questions will be the
    scalar.
    If \code{beta.start} is set to a matrix, the number
    of rows should be the number of questions and the number of columns
    should be 2.
    The first column
    will be the starting values for \eqn{\alpha_j} and the second column
    will be the starting values for \eqn{\beta_j}.
    Since the parameter values are constrained to be positive, the starting
    values should be also positive.
    The default is \code{1}.
  }

  \item{tau.start}{starting values for the cut points in the response
    model. If \code{NA}, the function generates the starting values
    so that each interval between the cut points is 0.5.
    If \code{tau.start} is set to a matrix, the number of rows should be
    the same as the number of questions and the number of columns should
    be the maximum value of the number of categories in the responses.
    The first cut point for each question should be set to 0 while the
    last one set to the previous cut point plus 1000.
    The default is \code{NA}.
  }

  \item{lambda.start}{starting values for the coefficients in the
    support parameter model, \eqn{\lambda_jk}{lambda_jk}.
    If \code{lambda.start} is set to a scalar, the starting values for
    all coefficients will be the scalar.
    If \code{lambda.start} is set to a matrix, the number of rows should
    be the number of the individual level covariates (plus the number of villages, if the model
    is hierarchical),
    and the number of columns should be the number of
    endorsers (times the number of questions, if the model is with varying lambdas).
    The default is \code{0}.
  }

  \item{omega2.start}{starting values for the variance of the support
    parameters, \eqn{\omega_{jk}^{2}}{omega_{jk}^2}.
    If \code{omega2.start} is set to a scalar, the starting values for
    \eqn{omega_{jk}^{2}}{omega_{jk}^2} will be the diagonal matrix with
    the diagonal elements set to the scalar.
    If \code{omega2.start} is set to a matrix, the number of rows should
    be the number of questions, 
    while the number of columns should be the same as the number of
    endorsers.
    The default is \code{.1}.
  }

  \item{theta.start}{starting values for the means of the
    \eqn{\lambda_{jk}}{lambda_jk} for each endorser.
    If \code{theta.start} is set to a scalar, the starting values for
    all parameters will be the scalar.
    If \code{theta.start} is set to a matrix, the number of rows should
    be the number of endorsers and the number of columns should be the
    dimension of covariates.
    The default is \code{0}.
  }

  \item{phi2.start}{starting values for the covariance matrices of the
    coefficients
    of the support parameters, \eqn{\Phi_{k}}{Phi_k}.
    \eqn{\Phi_{k}}{Phi_k} is assumed to be a diagonal matrix.
    If \code{phi2.start} is set to a scalar, the starting values for
    all covariance matrices will be the same diagonal matrix with the
    diagonal elements set to the scalar.
    If \code{phi2.start} is set to a vector, the length should be the
    number of endorsers times the dimension of covariates.
    The default is \code{.1}.
  }

  \item{kappa.start}{starting values for the coefficients on village level covariates in the
    support parameter model, \eqn{\kappa_k}{kappa_k}.
    If \code{kappa.start} is set to a scalar, the starting values for
    all coefficients will be the scalar.
    If \code{kappa.start} is set to a matrix, the number of rows should
    be the number of the village level covariates,
    and the number of columns should be the number of
    endorsers (times the number of questions, if the varying-lambda
    model is fitted).
    The default is \code{0}.
  }

  \item{psi2.start}{starting values for the variance of the village random intercepts in the support
    parameter model, \eqn{\psi_{k}^{2}}{psi_{k}^2}.
    If \code{psi2.start} is set to a scalar, the starting values for
    \eqn{\psi_{k}^{2}}{psi_{k}^2} will be the diagonal matrix with
    the diagonal elements set to the scalar.
    If \code{psi2.start} is set to a vector, its length should
    be the number of endorsers (times the number of questions, if the
    varying-lambda model is fitted).
    The default is \code{.1}.
  }

  \item{delta.start}{starting values for the coefficients on individual level covariates in the ideal
    point model. Will be used only if \code{covariates = TRUE}.
    If \code{delta.start} is set to a scalar, the starting values for
    all coefficients will be the scalar.
    If \code{delta.start} is set to a vector, the length should be the
    dimension of covariates.
    The default is \code{0}.
  }

  \item{zeta.start}{starting values for the coefficients on village level covariates in the ideal
    point model. Will be used only if \code{covariates = TRUE}.
    If \code{zeta.start} is set to a scalar, the starting values for
    all coefficients will be the scalar.
    If \code{zeta.start} is set to a vector, the length should be the
    dimension of covariates.
    The default is \code{0}.
  }

  \item{rho2.start}{numeric. starting values for the variance of the village random intercepts in the ideal point
    model, \eqn{\rho^{2}}{rho^2}. The default is \code{1}.
  }

  \item{mu.beta}{the mean of the independent Normal prior on the
    question related parameters. Can be either a scalar or a matrix of
    dimension the number of questions times 2.
    The default is \code{0}.
  }

  \item{mu.x}{the mean of the independent Normal prior on the
    question related parameters. Can be either a scalar or a vector of
    the same length as the number of observations.
    The default is \code{0}.
  }

  \item{mu.theta}{the mean of the independent Normal prior on the
    mean of the coefficients in the support parameter model.
    Can be either a scalar or a vector of
    the same length as the dimension of covariates.
    The default is \code{0}.
  }

  \item{mu.kappa}{the mean of the independent Normal prior on the
    coefficients of village level covariates. Can be either a scalar or a matrix of
    dimension the number of covariates times the number of endorsers.
    If auxiliary information is included, the value of \code{mu.kappa} 
    will be computed for each group such that the prior probability of
    the support parameter taking a positive value is equal to the known
    value of \code{h}.
    The default is \code{0}.
  }

  \item{mu.delta}{the mean of the independent Normal prior on the
    the coefficients in the ideal point model.
    Can be either a scalar or a vector of
    the same length as the dimension of covariates.
    The default is \code{0}.
  }

  \item{mu.zeta}{the mean of the independent Normal prior on the
    the coefficients of village level covariates in the ideal point model.
    Can be either a scalar or a vector of
    the same length as the dimension of covariates.
    The default is \code{0}.
  }

  \item{precision.beta}{the precisions (inverse variances) of the
    independent Normal prior on the
    question related parameters. Can be either a scalar or
    a 2 \eqn{\times}{by} 2 diagonal matrix.
    The default is \code{0.04}.
  }

  \item{precision.x}{scalar. The known precision of the
    independent Normal distribution on the
    ideal points.
    The default is \code{1}.
  }

  \item{precision.theta}{the precisions of the
    independent Normal prior on the means of the coefficients
    in the support parameter model. Can be either a scalar or
    a vector of the same length as the dimension of covariates.
    The default is \code{0.04}.
  }

  \item{precision.kappa}{the precisions of the
    independent Normal prior on the coefficients of village level covariates
    in the support parameter model. Can be either a scalar or
    a vector of the same length as the dimension of covariates.
    If auxiliary information is included, the value of \code{precision.kappa} 
    will be fixed to \code{100000}.
    The default is \code{0.04}.
  }

  \item{precision.delta}{the precisions of the
    independent Normal prior on the the coefficients
    in the ideal point model. Can be either a scalar or
    a square matrix of the same dimension as the dimension of
    covariates.
    The default is \code{0.04}.
  }

  \item{precision.zeta}{the precisions of the
    independent Normal prior on the the coefficients of village level covariates
    in the ideal point model. Can be either a scalar or
    a square matrix of the same dimension as the dimension of
    covariates.
    The default is \code{0.04}.
  }

  \item{s0.omega2}{scalar. The scale of the independent scaled
    inverse- chi-squared    
    prior for the variance parameter in the support parameter model.
    The default is \code{1}.
  }

  \item{nu0.omega2}{scalar. The degrees of freedom of the independent
    scaled inverse-chi-squared    
    prior for the variance parameter in the support parameter model.
    If auxiliary information is included, the value of \code{nu0.omega2} 
    will be fixed to the default.
    The default is \code{10}.
  }

  \item{s0.phi2}{scalar. The scale of the independent
    scaled inverse-chi-squared
    prior for the variances of the coefficients in
    the support parameter model.
    The default is \code{1}.
  }

  \item{nu0.phi2}{scalar. The degrees of freedom of the independent
    scaled
    inverse-chi-squared
    prior for the variances of the coefficients in
    the support parameter model.
    The default is \code{10}.
  }

  \item{s0.psi2}{scalar. The scale of the independent
    scaled inverse-chi-squared
    prior for the variances of the village random intercepts in
    the support parameter model.
    The default is \code{1}.
  }

  \item{nu0.psi2}{scalar. The degrees of freedom of the independent
    scaled
    inverse-chi-squared
    prior for the variances of the village random intercepts in
    the support parameter model.
    The default is \code{10}.
  }

  \item{s0.sig2}{scalar. The scale of the independent
    scaled inverse-chi-squared
    prior for the variance parameter in
    the ideal point model.
    The default is \code{1}.
  }

  \item{nu0.sig2}{scalar. The degrees of freedom of the independent
    scaled
    inverse-chi-squared
    prior for the variance parameter in the ideal point model.
    The default is \code{400}.
  }

  \item{s0.rho2}{scalar. The scale of the independent
    scaled inverse-chi-squared
    prior for the variances of the village random intercepts in
    the ideal point model.
    The default is \code{1}.
  }

  \item{nu0.rho2}{scalar. The degrees of freedom of the independent
    scaled
    inverse-chi-squared
    prior for the variances of the village random intercepts in
    the ideal point model.
    The default is \code{10}.
  }

  \item{MCMC}{the number of iterations for the sampler. The default is
    \code{20000}.
  }

  \item{burn}{the number of burn-in iterations for the sampler. The
    default is \code{1000}.
  }

  \item{thin}{the thinning interval used in the simulation. The default
    is \code{1}.
  }

  \item{mh}{logical. If \code{TRUE}, the Metropolis-Hastings algorithm
    is used to sample the cut points in the response model. The default is
    \code{TRUE}.
  }

  \item{prop}{a positive number or a vector consisting of positive
  numbers. The length of the vector should be the same as the number of
  questions. This argument sets proposal variance for the
  Metropolis-Hastings algorithm in sampling the cut points of the
  response model. The default is \code{0.001}.}

  \item{x.sd}{logical. If \code{TRUE}, the standard deviation of the
  ideal points in each draw will be stored. If \code{FALSE}, a sample of
  the ideal points will be stored. \emph{NOTE: Because storing a sample
  takes an enormous amount of memory, this option should be
  selected only if the chain is thinned heavily or the data have a small
  number of observations.}
  }

  \item{tau.out}{logical. A switch that determines whether or not to
    store the cut points in the response model. The default is
    \code{FALSE}.
  }
  
  \item{s.out}{logical. If \code{TRUE}, the support parameter for each
    respondent and each question will be stored.
    The default is \code{FALSE}. \emph{NOTE: Because storing a sample
    takes an enormous amount of memory, this option should be
    selected only if the chain is thinned heavily or the data have a small
    number of observations.}
  }

  \item{omega2.out}{logical. If \code{TRUE}, the variannce parameter of the support
    parameter model will be stored.
    The default is \code{TRUE}.
  }

  \item{phi2.out}{logical. If \code{TRUE}, the variannce parameter of the model
    for the coefficients in the support parameter model will be stored.
    The default is \code{TRUE}.
  }

  \item{psi2.out}{logical. If \code{TRUE}, the variance of the village random intercepts 
    in the support parameter model  will be stored.
    The default is \code{TRUE}.
  }

  \item{verbose}{logical. A switch that determines whether or not to
    print the progress of the chain and Metropolis acceptance ratios for
    the cut points of the response model. The default is
    \code{TRUE}.
  }

  \item{seed.store}{logical.  If \code{TRUE}, the seed will be stored in order
    to update the chain later.  The default is \code{FALSE}.
  }

  \item{update}{logical.  If \code{TURE}, the function is run to update a chain.
    The default is \code{FALSE}.
  }

  \item{update.start}{list.  If the function is run to update a chain, the output
    object of the previous run should be supplied.  The default is \code{NULL}.
  }
}

\details{
  The model takes the following form:
  
  Consider an endorsement experiment where we wish to measure the level
  of support for \eqn{K} political
  actors.  In the survey, respondents are asked whether or
  not they support each of \eqn{J} policies chosen by researchers.
  Let \eqn{Y_{ij}} represent respondent
  \eqn{i}'s answer to the survey question regarding policy \eqn{j}.
  Suppose that the response variable \eqn{Y_{ij}} is the ordered factor
  variable taking one of \eqn{L_{j}} levels, i.e., \eqn{Y_{ij} \in \{0,
    1, \dots,
  L_{j} - 1\}} where \eqn{L_{j} > 1}.  We assume that a greater value of
  \eqn{Y_{ij}}
  indicates a greater level of support for policy \eqn{j}.
  We denote an \eqn{M} dimensional vector of the observed
  characteristics of respondent \eqn{i} by \eqn{Z_i}.

  In the experiment, we
  randomly assign one of \eqn{K} political actors as an endorser to
  respondent \eqn{i}'s question regarding policy \eqn{j} and denote this
  treatment variable by \eqn{T_{ij} \in \{0,1,\dots,K\}}.  We use
  \eqn{T_{ij}=0} to represent the control observations where no
  political endorsement is attached to the question. Alternatively, one
  may use the endorsement by a neutral actor as the control group.
  
  The model for the response variable, \eqn{Y_{ij}}, is given by,
  \deqn{Y_{ij} = l \; {\rm if} \; \tau_{l} < Y^{*}_{ij} \le
    \tau_{l + 1},
  }{%
    Y_ij = l if \tau_l < Ystar_ij \le \tau_(l+1),%
  }
  \deqn{Y^{*}_{ij} \; | \; T_{ij} = k \sim \mathcal{N}(- \alpha_{j} +
    \beta_{j} (x_{i} + s_{ijk}), \; I)
  }{%
    Ystar_ij | T_ij = k ~ N(- \alpha_j + \beta_j (x_i + s_(ij, T_ij)), I)%
  }
  where \eqn{l \in \{0, 1, \dots, L_{j} \}, \tau_{0} = -\infty <
    \tau_{1} = 0 < \tau_{2} < \dots < \tau_{L_{j}} = \infty}.
  \eqn{\beta_j}'s are assumed to be positive.

  The model for the support parameter, \eqn{s_{ijk}}, is given by
  if \eqn{T_{ij} \neq 0},
  \deqn{
    s_{ijk} \sim \mathcal{N}(Z_i^{T} \lambda_{jk}, \; \omega_{jk}^2)
  }{%
    s_ijk ~ N(Z_i^T \lambda_jk, \omega_jk^2)%
  }
  with covariates, and
  \deqn{
    s_{ijk} \sim \mathcal{N}(\lambda_{jk}, \; \omega_{jk}^2),
  }{%
    s_ijk ~ N(\lambda_jk, \omega_jk^2)%
  }
  without covariates, for \eqn{j = 1, \dots, J, \; k = 1, \dots, K},
  and if \eqn{T_{ij} = 0, \; s_{ijk} = 0}{T_ij = 0, s_ijk = 0}.

  The \eqn{\lambda}'s in the support parameter model are modeled in the
  following hierarchical manner,
  \deqn{
    \lambda_{jk} \sim \mathcal{N}(\theta_k, \; \Phi_k)
  }{%
    \lambda_jk ~ N(\theta_k, \Phi_k)%
  }  
  for \eqn{k = 1, \dots, K}.


  If you set \code{identical.lambda = FALSE} and \code{hierarchical = TRUE},
  the model for \eqn{s_{ijk}} is if \eqn{T_{ij} \neq 0},
  \deqn{
    s_{ijk} \sim \mathcal{N}(\lambda^{0}_{jk, village[i]} + Z_i^{T} \lambda_{jk}, \; \omega_{jk}^2)
  }{%
    s_ijk ~ N(\lambda^{0}_{jk, village[i]} + Z_i^T \lambda_jk, \omega_jk^2)%
  }
  and
  \deqn{
    \lambda^{0}_{jk, village[i]} \sim \mathcal{N}(V_{village[i]}^{T} \kappa_{jk}, \; \psi_{jk}^2)
  }{%
    lambda^{0}_{jk, village[i]} ~ N( V_village[i]^T \kappa_jk, \psi_jk^2)%
  }
  for \eqn{k = 1, \dots, K} and \eqn{j = 1, \dots, J}.  In addition,
  \eqn{\lambda} and \eqn{\kappa} are modeled in the following
  hierarchical manner,
  \deqn{
    \lambda^{*}_{jk} \sim \mathcal{N}(\theta_k, \; \Phi_k)
  }{%
    \lambda*_jk ~ N(\theta_k, \Phi_k)%
  }  
  for \eqn{k = 1, \dots, K}, where \eqn{\lambda^{*}_{jk} =
  (\lambda^{T}_{jk}, \kappa^{T}_{jk})^{T}}.


  If you set \code{identical.lambda = TRUE} and \code{hierarchical = TRUE},
  the model for \eqn{s_{ijk}} is if \eqn{T_{ij} \neq 0},
  \deqn{
    s_{ijk} \sim \mathcal{N}(\lambda^{0}_{k, village[i]} + Z_i^{T} \lambda_{k}, \; \omega_{k}^2)
  }{%
    s_ijk ~ N(\lambda^{0}_{k, village[i]} + Z_i^T \lambda_k, \omega_k^2)%
  }
  and
  \deqn{
    \lambda^{0}_{k, village[i]} \sim \mathcal{N}(V_{village[i]}^{T} \kappa_{k}, \; \psi_{k}^2)
  }{%
    lambda^{0}_{k, village[i]} ~ N( V_village[i]^T \kappa_k, \psi_k^2)%
  }
  for \eqn{k = 1, \dots, K}.
  
  If the covariates are included in the model, the model for the
  ideal points is given by
  \deqn{
    x_{i} \sim \mathcal{N}(Z_{i}^{T} \delta, \; \sigma_{x}^{2})
  }{%
    x_i ~ N(Z_i^T \delta, \sigma_x^2)%
  }
  for \eqn{i = 1, \dots, N} where \eqn{\sigma_x^2} is a known prior
  variance.

  If you set \code{hierarchical = TRUE},
  the model is
  \deqn{
    x_{i} \sim \mathcal{N}(\delta^{0}_{village[i]} + Z_i^{T} \delta, \; \sigma^2)
  }{%
    x_i ~ N(\delta^{0}_village[i] + Z_i^T \delta, \sigma^2)%
  }
  and
  \deqn{
    \delta^{0}_{village[i]} \sim \mathcal{N}(V_{village[i]}^{T} \zeta, \; \rho^2)
  }{%
    delta^{0}_{village[i]} ~ N( V_village[i]^T \zeta, \rho^2)%
  }
  for \eqn{k = 1, \dots, K}.

  
  Finally, the following independent prior distributions are placed on
  unknown parameters,
  \deqn{
    \alpha_j \sim \mathcal{N}(\mu_\alpha, \; \sigma_\alpha^2)
  }{%
    \alpha_j ~ N(\mu_\alpha, \sigma_\alpha^2)%
  }
  for \eqn{j = 1, \dots, J},
  \deqn{
    \beta_j \sim \mathcal{TN}_{\beta_j > 0}(\mu_\beta, \; \sigma_\beta^2)
  }{%
    \beta_j ~ TN_{\beta_j > 0}(\mu_\beta, \sigma_\beta^2)%
  }
  for \eqn{j = 1, \dots, J},
  \deqn{
    \delta \sim \mathcal{N}(\mu_\delta, \; \Sigma_\delta),
  }{%
    \delta ~ N(\mu_\delta, \Sigma_\delta),%
  }
  \deqn{
    \theta_k \sim \mathcal{N}(\mu_\theta, \; \Sigma_\theta)
  }{%
    \theta_k \sim \mathcal{N}(\mu_\theta, \Sigma_\theta)%
  }
  for \eqn{k = 1, \dots, K},
  \deqn{
    \omega_{jk}^2 \sim  {\rm Inv-}\chi^{2}(\nu_{\omega}^0, \;
    s_{\omega}^0)
  }{%
    \omega_jk^2 ~ Inv-\chi^2(\nu_\omega^0, s_\omega^0)
  }
  for \eqn{j = 1, \dots, J} and \eqn{k = 1, \dots, K}, and
  \deqn{  
    {\rm diag}(\Phi_k) \sim  {\rm Inv-}\chi^{2}(\nu_{\Phi}^0, \;
    s_{\Phi}^0)
  }{%
    diag(\Phi_k) ~ Inv-\chi^2(\nu_\Phi^0, s_\Phi^0)
  }
  for \eqn{k = 1, \dots, K}, where \eqn{\Phi_k} is assumed to be a
  diagonal matrix.
}

\value{
  A list containing the following elements:

  \item{beta}{an mcmc object. A sample from the posterior distribution
    of \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
  }

  \item{x}{If \code{x.sd = TRUE}, a vector of the standard deviation of
    the ideal points in each draw. If \code{x.sd = FALSE}, an mcmc
    object that contains a sample from the posterior distribution of the
    ideal points.
  }

  \item{s}{If \code{s.out = TRUE}, an mcmc object that contains a sample
    from the posterior distribution of \eqn{s_{ijk}}{s_ijk}.  Variable
    names are \code{s(observation index)(question index)}.
  }

  \item{delta}{If \code{covariates = TRUE}, an mcmc object that contains
    a sample from the posterior distribution of \eqn{\delta}{delta}.
  }

  \item{tau}{If \code{tau.out = TRUE}, an mcmc object that contains a
    sample from the posterior distribution of \eqn{\tau}{tau}.
  }

  \item{lambda}{an mcmc object. A sample from the posterior distribution
    of \eqn{\lambda}{lambda}.  Variable names are:
    \code{lambda(question index)(group index).(covariate index)}.
  }

  \item{theta}{an mcmc object. A sample from the posterior distribution
    of \eqn{\theta}{theta}.
  }

  \item{kappa}{an mcmc object.
  }

  \item{zeta}{an mcmc object.
  }

  Note that the posterior sample of all parameters are NOT
  standardized. In making posterior inference, each parameter should be
  divided by the standard deviation of x (in the default setting, it is
  given as "x") or by \eqn{\sigma^2}{sigma} (in the default setting, it
  is given as "sigma2").

  Also note that \eqn{\alpha}{alpha} and the intercept in
  \eqn{\delta}{delta} (or, if the model is hierarchical, the intercept
  in \eqn{\zeta}{zeta}) are not identified.  Instead,
  \deqn{  
    - \alpha + \beta * \delta_0
  }{%
    - alpha + beta * delta_0
  }
  or, if the model is hierarchical,
  \deqn{  
    - \alpha + \beta * \zeta_0
  }{%
    - alpha + beta * zeta_0
  }
  is identified after either of the above standardization, where
  \eqn{\delta_0}{delta_0} and \eqn{\zeta_0}{zeta_0} denote the
  intercepts.
  
  When using the auxiliary data functionality, the following objects are included:

  \item{aux}{logical value indicating whether estimation incorporates auxiliary moments
  }
  
  \item{nh}{integer count of the number of auxiliary moments
  }
}

\references{
  Bullock, Will, Kosuke Imai, and Jacob N. Shapiro. (2011)
  \dQuote{Statistical Analysis of Endorsement Experiments: Measuring
  Support for Militant Groups in Pakistan,} \emph{Political Analysis},
  Vol. 19, No. 4 (Autumn), pp.363-384.
}

\examples{
\dontrun{
data(pakistan)

Y <- list(Q1 = c("Polio.a", "Polio.b", "Polio.c", "Polio.d", "Polio.e"),
          Q2 = c("FCR.a", "FCR.b", "FCR.c", "FCR.d", "FCR.e"),
          Q3 = c("Durand.a", "Durand.b", "Durand.c", "Durand.d",
                 "Durand.e"),
          Q4 = c("Curriculum.a", "Curriculum.b", "Curriculum.c",
                 "Curriculum.d", "Curriculum.e"))



## Varying-lambda non-hierarchical model without covariates
endorse.out <- endorse(Y = Y, data = pakistan, identical.lambda = FALSE,
                       covariates = FALSE, hierarchical = FALSE)

## Varying-lambda non-hierarchical model with covariates
indiv.covariates <- formula( ~ female + rural)
endorse.out <- endorse(Y = Y, data = pakistan, identical.lambda = FALSE,
                       covariates = TRUE, formula.indiv = indiv.covariates,
                       hierarchical = FALSE)

## Common-lambda non-hierarchical model with covariates
indiv.covariates <- formula( ~ female + rural)
endorse.out <- endorse(Y = Y, data = pakistan, identical.lambda = TRUE,
                       covariates = TRUE, formula.indiv = indiv.covariates,
                       hierarchical = FALSE)

## Varying-lambda hierarchical model without covariates
div.data <- data.frame(division = sort(unique(pakistan$division)))
div.formula <- formula(~ 1)
endorse.out <- endorse(Y = Y, data = pakistan, data.village = div.data,
                       village = "division", identical.lambda = FALSE,
                       covariates = FALSE, hierarchical = TRUE,
                       formula.village = div.formula)

## Varying-lambda hierarchical model with covariates
endorse.out <- endorse(Y = Y, data = pakistan, data.village = div.data,
                       village = "division", identical.lambda = FALSE,
                       covariates = TRUE,
                       formula.indiv = indiv.covariates,
                       hierarchical = TRUE,
                       formula.village = div.formula)

## Common-lambda hierarchical model without covariates
endorse.out <- endorse(Y = Y, data = pakistan, data.village = div.data,
                       village = "division", identical.lambda = TRUE,
                       covariates = FALSE, hierarchical = TRUE,
                       formula.village = div.formula)

## Common-lambda hierarchical model with covariates
endorse.out <- endorse(Y = Y, data = pakistan, data.village = div.data,
                       village = "division", identical.lambda = TRUE,
                       covariates = TRUE,
                       formula.indiv = indiv.covariates,
                       hierarchical = TRUE,
                       formula.village = div.formula)
}
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu}

  Yuki Shiraito, Department of Politics, Princeton University
  \email{shiraito@Princeton.Edu}.
}

