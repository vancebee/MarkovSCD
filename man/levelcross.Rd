\name{levelcross}

\alias{levelcross}

\title{levelcross}

\description{Calculate an empirical level-crossing function to assist in the selection of state boundaries.}

\usage{levelcross(tseries, evalpts = NULL, npts = 100, lag=1)}


\arguments{
  \item{tseries}{A vector of values defining a time series.  The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{evalpts}{Optional; the points at which the level-crossing function is evaluated can be supplied.  If not, \code{evalpts} will be automatically determined according to the min/max of \code{tseries} and \code{npts}.}

  \item{npts}{If \code{evalpts} are not supplied, \code{npts} points evenly-spaced between the minimum and maximum of \code{tseries} will be used.}

  \item{lag}{Number of observations separating points that will be assessed for level crossing. }

}

\value{
  A list consisting of two elements:

  \item{evalpts}{Either the supplied or calculated points at which the level-crossing function is evaluated.}

  \item{lvlcrs}{A vector where each element corresponds to the proportion of the data that crossed the corresponding element of \code{evalpts}.}
}

