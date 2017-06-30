\name{validitycheck}

\alias{validitycheck}

\title{validitycheck}

\description{Calculate the Frobenious norm of the delta matrix for variosu stratification points}

\usage{validitycheck(tseries1, tseries2, ilength1, ilength2,
statebounds,lag, nruns = 10000, CI = 0.95, mineffect = 0.05)}


\arguments{

    \item{tseries1,
  tseries2}{Vectors of values defining two time series. The last observation in \code{tseries1} should be considered to immediately precede the first observation of \code{tseries2}. The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{ilength1,
  ilength2}{The number of observations required for the dynamics in \code{tseries1,tseries2} to converge. Although not required, these values will likely come from calling \code{dynamicsconv}.}

  \item{statebounds}{Boundaries delineating the states to be considered.}

  \item{lag}{The number of measurements separating each observation and the lagged observation. Transitions are counted between the ith and (i+\code{lag})th measurements.}

  \item{nruns}{Number of bootstrap runs to determine significance.}

  \item{CI}{Width of confidence interval to be used to determine significance.  Must be within (0,1).}

  \item{mineffect}{Absolute value of smallest effect size that will be considered significant.}
}

\value{
  A list consisting of two elements:

  \item{norm}{A vector where each element is the norm of the delta matrix associated with an iteration.}

  \item{diagconfig}{A matrix where the top, middle, and bottom rows correspond with the mean of the superdiagonal, diagonal, and subdiagonal entries of an iteration associated with each column.}

}

