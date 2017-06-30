lageval = function(tseries,statebounds,lagrange){##Function to find transition prob matrix

  ###########################Check formatting######################
  if(!is.numeric(tseries)) stop("Time-series is not a numeric vector")

  if(!is.numeric(statebounds)) stop("statebounds is not a numeric vector")

  if(length(statebounds) <= 2) stop("Must have at least 3 statebounds (2 bins)")

  ######Calc transition matrix for each lag being considered######
  nlags = length(lagrange)

  tmlist = list()
  for(ii in 1:nlags){
    tm =  transmat(tseries,statebounds,lagrange[ii])
    tmlist[[ii]] = tm$prob
  }

  #For each state, extract the diagonal element associated with####
  #each of the transition matrices produced in previous step#######
  #################################################################
  nstates = length(statebounds)-1
  nm = NULL
  for(ii in 1:nstates) nm = c(nm,paste("S",ii,sep=""))

  diagbylag = list()
  for(ii in 1:nstates){
    diagbylag[[ii]] = unlist(lapply(tmlist, function(z) z[ii,ii]))
  }
  names(diagbylag) = nm

  ###Return lag range and list diag probs##########################
  rtn = list(lagrange,diagbylag)
  names(rtn) = c("lagrange", "diagbylag")

  return(rtn)
}
