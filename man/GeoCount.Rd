\name{GeoCount}
\alias{GeoCount}
\title{
  Counting Violent Incidents around Villages
}
\description{
  This function calculates the number of violent incidents within a specified
  distance around specified points (villages).
}
\usage{
GeoCount(x, y, distance, x.latitude = "latitude",
                      x.longitude = "longitude", y.latitude = "latitude",
                      y.longitude = "longitude")
}
\arguments{

 \item{x}{data frame containing the longitude and the latitude of villages.
 }

 \item{y}{data frame containing the longitude and the latitude of violent
   incidents.
 }

 \item{distance}{numeric. The distance from villages in kilometers.
 }

 \item{x.latitude}{character.  The variable name for the latitude in \code{x}.
 }

 \item{x.longitude}{character.  The variable name for the longitude in \code{x}.
 }

 \item{y.latitude}{character.  The variable name for the latitude in \code{y}.
 }

 \item{y.longitude}{character.  The variable name for the longitude in \code{y}.
 }

}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu}

  Yuki Shiraito, Department of Politics, Princeton University
  \email{shiraito@Princeton.Edu}.
}

