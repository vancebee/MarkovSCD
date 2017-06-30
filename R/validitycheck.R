validitycheck = function(tseries1, tseries2, ilength1, ilength2, statebounds,lag, nruns = 10000, CI = 0.95, mineffect = 0.05){

  ###########################Check formatting######################
  if(!is.numeric(tseries1)) stop("Time-series 1 is not a numeric vector")

  if(!is.numeric(tseries2)) stop("Time-series 2 is not a numeric vector")

  #################################################################
  #Combine two time series into 1 and calc indices defining
  #intervals to be considered
  #################################################################
  nn1 = length(tseries1)
  nn2 = length(tseries2)

  tseries = c(tseries1,tseries2)
  nn = length(tseries)
  ilength = ilength1+ilength2

  #Select window so there is a matchup witth BL/TX Interface after 4 iterations
  delta = round((nn1-ilength1)/4)
  #Caclu number of iterations before reaching the end
  nits=1
  while(ilength + nits*delta < nn) nits = nits+1

  #################################################################
  #Find norm of each delta matrix
  #################################################################
  storedelta = list()
  for(ii in 1:nits){
    ts1l = (ii-1)*delta+1; ts1u = ts1l+ilength1
    ts2l = ts1u+1        ; ts2u = ts2l+ilength2

    ts1 = tseries[ts1l:ts1u]
    ts2 = tseries[ts2l:ts2u]

    tm1 = transmat(ts1,statebounds,lag)
    tm2 = transmat(ts2,statebounds,lag)

    dt = deltatrans(tm1,tm2,nruns,CI,mineffect)

    storedelta[[ii]] = dt$prettydelta
  }

  fnorm = unlist(lapply(storedelta,function(x) norm(x,type = "f")))

  #################################################################
  #Find the average super/diag/sub for each delta matrix
  #################################################################
  super = apply(sapply(storedelta,function(z) diag(z[-nrow(z),-1])),2,mean)
  dgl = apply(sapply(storedelta,diag),2,mean)
  sub = apply(sapply(storedelta,function(z) diag(z[-1,-ncol(z)])),2,mean)

  MM = matrix(c(super,dgl,sub),nrow = 3, byrow = T)
  rownames(MM) = c("super","diag","sub")

  #################################################################
  #Return norm of differences and range of indices
  #################################################################
  rtn = list(fnorm,MM)
  names(rtn) = c("norm","diagconfig")

  return(rtn)
}
