\name{sensitivity}

\alias{sensitivity}

\title{sensitivity}

\description{Produces delta matrices comparing two different phases for various combinations of lag and state boundaries.}

\usage{sensitivity(tseries1, tseries2, stbdyrange,lagrange)}


\arguments{
  \item{tseries1}{A vector of values defining a time series.  The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{tseries2}{A vector of values defining a time series.  The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{stbdyrange}{A list of vectors. Each list element provides boundaries for delineating the states.}

  \item{lagrange}{A vector indicating which lags, in units of number of observations, are to be considered.}
}

\value{
  A list consisting of three elements:

  \item{stbdyrange}{The \code{lagrange} that was input to the function. For reference in \code{deltamats}, each element of the list has been named "Bdy1", "Bdy2",...}

  \item{lagrange}{The \code{lagrange} that was input to the function.}

  \item{deltamats}{A list of length \code{length(stbdyrange)*length(lagrange)}.  Each element consists of a \code{prettydelta} matrix summarizing the difference between the dynamics of \code{tm1} and \code{tm2}.  For reference, each element is named "BdyX_by_LagX", where X comes from the names supplied in \code{stbdyrange} and Y is a value from \code{lagrange}.}
}

