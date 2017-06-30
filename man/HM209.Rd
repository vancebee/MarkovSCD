\name{HM209}

\alias{HM209}

\docType{data}

\title{Project Fresh Air Data from Home 209}

\description{Deindetified data from Home 209 of Project Fresh Air.  The dataframe consists of two variables:

1. a time series of air particle mass averages, each separated by 10 sec, and

2. an indicator variable indicating whether the observation is from the Baseline or Treatment phase of the intervention

}
\usage{data("HM209")}
\format{
  A data frame with 424519 observations on the following 2 variables.
  \describe{
    \item{\code{MassAve}}{a numeric vector}
    \item{\code{Phase}}{a character vector}
  }
}

\references{
For a full description of the trial, see:
Hughes S., J.  Bellettiere, B. Nguyen, S. Liles, N. Klepeis, P. Quintana, V. Berardi; S. Obayashi, S. Bradley, C. Hofstetter, M. Hovell. “Randomized Trial to Reduce Air Particle Levels in Homes of Smokers and Children.” PsyArXiv [Preprint].  Available from http://doi.org/10.17605/OSF.IO/M3V2Z
}
\examples{
data(HM209)
}
\keyword{datasets}
