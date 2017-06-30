\name{transmat}

\alias{transmat}

\title{transmat}

\description{Create an empirical transition matrix from time-series data. States are defined by different ranges of time series values.  The state of each observation (source) and state of the observation located some time after this observation (destination) are recorded and the number of transitions from each source to each destination state are calculated.}

\usage{transmat(tseries,statebounds,lag)}


\arguments{
  \item{tseries}{A vector of values defining a time series.  The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the \code{lag}.}

  \item{statebounds}{Boundaries delineating the states to be considered.}

  \item{lag}{The number of measurements separating each observation and the lagged observation. Transitions are counted between the ith and (i+\code{lag})th measurements.}
}

\value{
  A list consisting of two elements:

  \item{prob}{A transition matrix with source states as rows and destination states as columns. The (i,j)th entry corresponds to the empirical probability of transitioning from state i to state j.}

  \item{raw}{Similar to \code{prob} but elements are the raw number of transitions between states.}
}

