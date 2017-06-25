\name{predict.endorse} 
\alias{predict.endorse} 
\title{Predict Method for the Measurement Model of Political Support}
\description{ 
Function to calculate predictions and uncertainties of predictions from estimates from multivariate regression analysis of survey data with the item count technique.
} 
\usage{

  \method{predict}{endorse}(object, newdata, type = c("linear.s", "prob.support"), standardize = TRUE, ...) 

} 

\arguments{ 
\item{object}{Object of class inheriting from "endorse"}
\item{newdata}{An optional data frame containing data that will be used to make predictions from. If omitted, the data used to fit the regression are used.}
\item{type}{}
\item{standardize}{}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{\code{predict.endorse} produces predicted values, obtained by evaluating the regression function in the frame newdata (which defaults to \code{model.frame(object)}. If the logical \code{se.fit} is \code{TRUE}, standard errors of the predictions are calculated. Setting \code{interval} specifies computation of confidence intervals at the specified level or no intervals.

}


\value{
}

\seealso{\code{\link{endorse}} for model fitting}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu}

  Yuki Shiraito, Department of Politics, Princeton University
  \email{shiraito@Princeton.Edu}.
}

