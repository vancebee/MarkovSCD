\name{deltatrans}

\alias{deltatrans}

\title{deltatrans}

\description{Calculate the difference between two transition matrices and determine which differences are significant.}

\usage{deltatrans(tm1, tm2, nruns = 10000, CI = 0.95, mineffect = 0.05)}


\arguments{
  \item{tm1}{First transition matrix. Must be in format that is produced by a call to \code{transmat}}

  \item{tm2}{First transition matrix. Must be in format that is produced by a call to \code{transmat}}

  \item{nruns}{Number of bootstrap runs.}

  \item{CI}{Width of confidence interval to be used to determine significance.  Must be within (0,1).}

  \item{mineffect}{Absolute value of smallest effect size that will be considered significant.}
}

\value{
  A list consisting of three elements:

  \item{delta}{Matrix summarizing \code{tm2$prob}-\code{tm1$prob}}.

  \item{sigmat}{Matrix with a 1 for significant differences and 0 otherwise.}

  \item{prettydelta}{\code{delta} with only signficant values, rounded to two decimal places, printed.}
}

