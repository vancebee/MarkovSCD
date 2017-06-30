\name{mfpt}

\alias{mfpt}

\title{mfpt}

\description{Given a transition matrix, calculate mean time to reach one state for the first time from all other states.}

\usage{mfpt(tm, deststate)}


\arguments{
  \item{tm}{A transition probability matrix.}

  \item{deststate}{Index of the destination state.  The \code{deststate} row and column are removed from \code{tm} to calcualte mfpt.}
}

\value{
  A vector with one element for all states other than \code{deststate}.  The value represents the mfpt to reach \code{deststate} from this state.  A names attribute provides the index of states associated with the elements of the vector.
}

