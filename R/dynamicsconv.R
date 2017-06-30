dynamicsconv = function(tseries1, tseries2, nitvl, statebounds,lag, nruns = 10000, CI = 0.95, mineffect = 0.05){##Function to find transition prob matrix

  ###########################Check formatting######################
  if(!is.numeric(tseries1)) stop("Time-series 1 is not a numeric vector")

  if(!is.numeric(tseries2)) stop("Time-series 2 is not a numeric vector")

  #################################################################
  #Calc indices defining intervals to be considered
  #################################################################
  nn1 = length(tseries1)
  nn2 = length(tseries2)

  rng1 = round(seq(nn1-nn1/nitvl,1,-nn1/nitvl))
  if(rng1[length(rng1)]!=nn1) rng1[length(rng1)+1] = 1
  rng2 = round(seq(nn2/nitvl,nn2,nn2/nitvl))
  if(rng2[length(rng2)]!=nn2) rng2[length(rng2)+1] = nn2

  #################################################################
  #For each iteration, calcualte delta matrix
  #################################################################
  storedelta = list()
  for(ii in 1:length(rng1)){
    ts1 = tseries1[rng1[ii]:nn1]
    ts2 = tseries2[1:rng2[ii]]

    tm1 = transmat(ts1,statebounds,lag)
    tm2 = transmat(ts2,statebounds,lag)

    dt = deltatrans(tm1,tm2,nruns,CI,mineffect)

    storedelta[[ii]] = dt$prettydelta
  }

  #################################################################
  #Find norm of the difference between each matrix iteration and
  #matrix from using all data
  #################################################################
  fullmat = storedelta[[10]]
  fnorm = unlist(lapply(storedelta,function(x) norm(x-fullmat,type = "f")))

  #################################################################
  #Return norm of differences and range of indices
  #################################################################
  l1 = -cumsum(diff(c(nn1,rng1)))
  l2 = cumsum(diff(c(1,rng2)))

  rtn = list(fnorm,rng1,rng2,l1,l2)
  names(rtn) = c("normdiff","itvl1","itvl2","ilength1","ilength2")

  return(rtn)
}
