\name{lageval}

\alias{lageval}

\title{lageval}

\description{Evaluate the effects of different lags on the transition matrix.  For each lag considered, an empirical transition matrix is calculated and the diagonal elements from this matrix are extracted so that autocorrelation effects can be considered.}

\usage{lageval(tseries,statebounds,lagrange)}


\arguments{
  \item{tseries}{A vector of values defining a time series.  The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{statebounds}{Boundaries delineating the states to be considered.}

  \item{lagrange}{A vector indicating which lags, in units of number of observations, are to be considered.}
}

\value{
  A list consisting of two elements:

  \item{lagrange}{The \code{lagrange} that was input to the function.}

  \item{diagbylag}{A list of vectors. Each element of this list corresponds to the diagonal
  elements of one state (defined by \code{statebounds}). The vector that comprises each     element is the same length as \code{lagrange} and describes the diagonal element of the associated state at each lag under consideration.}
}

