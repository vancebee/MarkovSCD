\name{dynamicsconv}

\alias{dynamicsconv}

\title{dynamicsconv}

\description{Calcualte delta matrix for increasingly larger subsets of the data and compere to the delta matrix calcualted using all data via the Frobenious norm.}

\usage{dynamicsconv(tseries1, tseries2, nitvl, statebounds,lag,
nruns = 10000, CI = 0.95, mineffect = 0.05)}


\arguments{

  \item{tseries1,
  tseries2}{Vectors of values defining two time series. The last observation in \code{tseries1} should be considered to immediately precede the first observation of \code{tseries2}. The times at which these values were collected are not recorded, but knowledge of this feature will inform the selection of the lag.}

  \item{nitvl}{Number of evenly spaced intervals to stratify \code{tseries1} and \code{tseries2} into.  The first interval uses the
  last \code{length(tseries1)/nitvl})th of \code{tseries1} and the first \code{length(tseries2)/nitvl})th of \code{tseries2}. The second interval uses the
  last \code{2*length(tseries1)/nitvl})th of \code{tseries1} and the first \code{2*length(tseries2)/nitvl})th of \code{tseries2}, and so on.}

  \item{nruns}{Number of bootstrap runs to determine significance.}

  \item{CI}{Width of confidence interval to be used to determine significance.  Must be within (0,1).}

  \item{mineffect}{Absolute value of smallest effect size that will be considered significant.}
}

\value{
  A list consisting of five elements:

    \item{normdiff}{A vector with each element consisting of the norm of the difference in delta matrices between using all of the data and each interval under consideration.}

  \item{itvl1}{Vector of the indices of the intervals used for \code{tseries1}. For the ith iteration, indices \code{(itvl1[i]):length(tseries1)} were used to form a transition matrix.}

  \item{itvl2}{Vector of the indices of the intervals used for \code{tseries2}. For the ith iteration, indices \code{(1:itvl2[i])} were used to form a transition matrix.}

  \item{ilength1}{Length of each interval used for \code{tseries1}. It is the \code{cumsum} of the differences of \code{itvl1}.}

  \item{ilength2}{Length of each interval used for \code{tseries2}. It is the \code{cumsum} of the differences of \code{itvl2}.}
}

