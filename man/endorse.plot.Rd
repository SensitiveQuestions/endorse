\name{endorse.plot}
\alias{endorse.plot}
\title{
  Descriptive Plot of Endorsement Experiment Data
}
\description{
  This function creates a descriptive plot for a question in an
  endorsement experiment.
}
\usage{
endorse.plot(Y, data, scale, dk = 98, ra = 99, yaxis = NULL,
             col.seq = NA)
}
\arguments{

  \item{Y}{a character vector. List of the variable names for the
    responses to a question. Each variable name corresponds to each
    treatment status.
  }

  \item{data}{data frame containing the variables.
  }

  \item{scale}{an integer. The scale of the responses. The function
    assumes that the responses are coded so that \code{1} indicates the
    lowest support while the integer specified in this argument
    represents the highest support.
  }
  
  \item{dk}{an integer indicating the value of the response variable
    that is to be interpreted as ``Don't Know.'' Default is \code{98}.
  }

  \item{ra}{an integer indicating the value of the response variable
    that is to be interpreted as ``Refused.'' Default is \code{99}.
  }

  \item{yaxis}{a character vector of the same length as \code{Y}. The
    argument will be used for the label of the horizontal axis. The
    order should be the same as \code{Y}.
  }

  \item{col.seq}{a vector of colors for the bars or bar components. By
  default, a gradation of gray where the darkest indicates the highest
  support level.
  }

}

\value{
  A descriptive plot for the responses to a question.
}

\examples{

data(pakistan)

Y <- c("Polio.a", "Polio.b", "Polio.c", "Polio.d", "Polio.e")
yaxis <- c("Control", "Kashmir", "Afghan", "Al-Qaida", "Tanzeems")

endorse.plot(Y = Y, data = pakistan, scale = 5)

}

\author{
  Yuki Shiraito, Department of Political Science, University of Michigan
  \email{shiraito@umich.edu}.

  Kosuke Imai, Department of Government and Statistics, Harvard University
  \email{Imai@Harvard.Edu}, \url{https://imai.fas.harvard.edu/}
}

