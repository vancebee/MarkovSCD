transmat = function(tseries,statebounds,lag){##Function to find transition prob matrix

  ###############Check formatting#################################
  if(!is.numeric(tseries)) stop("Time-series is not a numeric vector")

  if(!is.numeric(statebounds)) stop("statebounds is not a numeric vector")

  if(!is.numeric(lag)) stop("Lag is not numeric")

  if(length(statebounds) <= 2) stop("Must have at least 3 statebounds (2 bins)")

  ##############Create lagged variable##########################
  nbds = length(statebounds)
  nn = length(tseries)
  ll = ifelse((1:nn)+lag<=nn, tseries[(1:nn)+lag], NA)

  #######Determine states of current and lagged variable########
  sourcebin = findInterval(tseries,statebounds)
  destbin = findInterval(ll,statebounds)

  #########Create raw and prob transition matrix###############
  nstates = nbds-1
  tm.raw =  matrix(0, nrow = nstates, ncol = nstates)

  for(ii in 1:nstates){
    for(jj in 1:nstates){
      tm.raw[ii,jj] = length(which(sourcebin == ii & destbin == jj))
    }
  }

  rs = rowSums(tm.raw)
  tm.prob = tm.raw/rs

  ######Name rows and columns of transition matrices#############
  nm = NULL
  for(ii in 1:nstates) nm = c(nm,paste("S",ii,sep=""))

  colnames(tm.prob) = nm; rownames(tm.prob) = nm
  colnames(tm.raw) = nm; rownames(tm.raw) = nm


  ##############Return list of matrices##########################
  rtn = list(tm.prob,tm.raw)
  names(rtn) = c("tm.prob", "tm.raw")

  return(rtn)
}
