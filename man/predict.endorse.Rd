\name{predict.endorse} 
\alias{predict.endorse} 
\title{Predict Method for the Measurement Model of Political Support}
\description{ 
Function to calculate predictions from a measurement model fitted to an
endorsement experiment data.
} 

\usage{\method{predict}{endorse}(object, newdata, type = c("prob.support",
"linear.s"), standardize = TRUE, ...)
} 

\arguments{ 
\item{object}{a fitted object of class inheriting from \code{"endorse"}}
\item{newdata}{an optional data frame containing data that will be used
  to make predictions from. If omitted, the data used to fit the
  regression are used.}
\item{type}{the type of prediction required.  The default is on the
  scale of the predicted probability of positive support; the
  alternative \code{"linear.s"} is on the scale of \eqn{s_{ijk}}.}
\item{standardize}{logical switch indicating if the predicted values on
  the scale of \eqn{s_{ijk}} are standardized so that its variance is
  one.}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{
  \code{predict.endorse} produces predicted support for political actors
  from a fitted \code{"endorse"} object.  If \code{newdata} is omitted
  the predictions are based on the date used for the fit.  Setting
  \code{type} specifies the type of predictions.  The default is
  \code{"prob.support"}, in which case the function computes the average
  predicted probability of positive support:
  \deqn{
    P(s_{ijk} > 0 \mid Z_i, \; \lambda_{j}, \; \omega_{j}) =
    \Phi \left( \frac{ Z_i^{T} \lambda_{j} }{ \omega_{j} } \right)
  }{%
    P(s_ijk > 0 | Z_i, \lambda_j, \omega_j) =
    \Phi(  Z_i^T \lambda_j / \omega_j )%
  }
  for each political group \eqn{k}. If \code{type} is set to be
  \code{"linear.s"}, the output is the predicted mean of support
  parameters:
  \deqn{
    E(s_{ijk} \mid Z_i, \; \lambda_{j}) = Z_i^{T} \lambda_{j}.
  }{%
    E(s_ijk | Z_i, \lambda_j) = Z_i^T \lambda_j.%
  }
  If the logical \code{standardize} is \code{TRUE}, the predicted mean
  of support is standardized by dividing by \eqn{\omega_j}.
}

\value{
  A \code{"mcmc"} object for predicted values.
}

\seealso{\code{\link{endorse}} for model fitting}

\author{
  Yuki Shiraito, Department of Politics, Princeton University
  \email{shiraito@Princeton.Edu}.
}

