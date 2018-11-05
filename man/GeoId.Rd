\name{GeoId}
\alias{GeoId}
\title{
  Getting Indices of Incidents around a specified point
}
\description{
  This function obtains the indices of incidents within a specified
  distance around a specified point.
}
\usage{
GeoId(x, y, distance, x.latitude = "latitude",
      x.longitude = "longitude", y.latitude = "latitude",
      y.longitude = "longitude")

}
\arguments{

 \item{x}{data frame containing the longitude and the latitude of a
   point.
 }

 \item{y}{data frame containing the longitude and the latitude of
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

\value{
  A vector containing the indices of \code{y} that are within
  \code{distance} kilometers around the point specified by \code{x}.  If
  there are multiple observations in \code{x}, the first row is used as
  the point.
}

\author{
  Yuki Shiraito, Department of Political Science, University of Michigan
  \email{shiraito@umich.edu}.
}

