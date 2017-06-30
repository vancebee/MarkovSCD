#Function to calculate and determine the significance between two transition matrices

levelcross = function(tseries, evalpts = NULL, npts = 100, lag=1){

  ###############Check formatting#################################
  if(!is.numeric(tseries)) stop("Time-series is not a numeric vector")

  #################################################################
  #If evalpts are not provided, create them based in min/max of data
  if(is.null(evalpts)){
    dr = diff(range(tseries))
    evalpts = seq(min(tseries)+dr/(npts+1),max(tseries)-dr/(npts+1),length.out = npts)
  }

  lvlcrs = NULL
  for(ii in 1:length(evalpts)){
    lvl = evalpts[ii]
    ltlvl = ifelse(tseries<lvl,1,0)

    aa = ltlvl[1:(length(ltlvl)-lag)]
    bb = ltlvl[(lag+1):length(ltlvl)]
    lvlcrs[ii] = length(which((aa-bb)!=0))
  }
  lvlcrs = lvlcrs/length(tseries)

  rtn = list(evalpts,lvlcrs)
  names(rtn) = c("evalpts","lvlcrs")

  return(rtn)
}



