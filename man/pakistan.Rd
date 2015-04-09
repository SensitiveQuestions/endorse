\name{pakistan}
\docType{data}
\alias{pakistan}
\title{Pakistan Survey Experiment on Support for Militant Groups}
\description{
  This data set is a subset of the data from the endorsement experiment
  conducted in Pakistan to study support for militant groups. The survey
  was implemented by Fair et al. (2009). It is also used by Bullock et
  al. (2011).
}
\usage{data(pakistan)}
\format{A data frame containing 5212 observations.
  The variables are:
  \itemize{
    \item \code{division}: division number.
    \item \code{edu}: education. 1 if ``illiterate''; 2 if
    ``primary''; 3 if ``middle''; 4 if ``matric''; 5 if ``intermediate
    (f.a/f.sc),'' ``graduate (b.a/b.sc.),'' or ``professionals (m.a /or
    other professional degree).''
    \item \code{inc}: approximate monthly income. 1 if less than 3000
    rupees; 2 if 3000 to 10,000 rupees; 3 if 10,001 to 15,000 rupees; 4
    if more than 15,000 rupees.
    \item \code{female}: 0 if male; 1 if female
    \item \code{rural}: 0 if rural; 1 if urban
    \item \code{Polio.a-e}: support for World Health Organization's plan of
    universal polio vaccinations in Pakistan. 5 indicates the highest
    support while 1 indicates the lowest support.
    \item \code{FCR.a-e}: support for the reform of the Frontier Crimes
    Regulation (FCR) governing the tribal areas. 5 indicates the highest
    support while 1 indicates the lowest support.
    \item \code{Durand.a-e}: support for using peace jirgas to resolve disputes
    over the Afghan border, the Durand Line. 5 indicates the highest
    support while 1 indicates the lowest support.
    \item \code{Curriculum.a-e}: support for the Government of Pakistan's plan
    of curriculum reforms in religious schools or \emph{madaris}. 5
    indicates the highest
    support while 1 indicates the lowest support.
  }
  For the response variables, endorsers are:
  \itemize{
    \item \code{varname.a}: control (no endorsement).
    \item \code{varname.b}: Pakistani militant groups in Kashmir.
    \item \code{varname.c}: Militants fighting in Afghanistan.
    \item \code{varname.d}: Al-Qaida.
    \item \code{varname.e}: Firqavarana Tanzeems.
  }
  
}
\source{Bullock, Will, Kosuke Imai, and Jacob
  N. Shapiro. 2011. Replication data
  for: Statistical analysis of endorsement experiments: Measuring
  support for militant groups in Pakistan. hdl:1902.1/14840. The
  Dataverse Network.
}
\references{
  Bullock, Will, Kosuke Imai, and Jacob N. Shapiro. (2011)
  \dQuote{Statistical Analysis of Endorsement Experiments: Measuring
  Support for Militant Groups in Pakistan,} \emph{Political Analysis},
  Vol. 19, No. 4 (Autumn), pp.363-384.

  Fair, Christin C., Neil Malhotra, and Jacob N. Shapiro. (2009)
  \dQuote{The Roots of Militancy: Explaining Support for Political
    Violence in Pakistan,} Working Paper, Princeton University.
}
\keyword{dataset}