sensitivity = function(tseries1, tseries2, stbdyrange,lagrange){##Function to find transition prob matrix

  ###########################Check formatting######################
  if(!is.numeric(tseries1)) stop("Time-series 1 is not a numeric vector")

  if(!is.numeric(tseries2)) stop("Time-series 2 is not a numeric vector")

  if(!class(stbdyrange) == "list") stop("A list of state boundaries must be provided")

  #################################################################
  ######Calc # of lags and # of boundary ranges to be considered ##
  #Name boundary Ranges############################################
  nlags = length(lagrange)

  nst = length(stbdyrange)
  for(ii in 1:nst) names(stbdyrange)[ii] = paste("Bdy",ii,sep="")

  #################################################################
  #For each combo of lag/state boundary, calc the difference between
  #two matrices and return the pretty delta matrix
  #Name boundary Ranges############################################
  ctr = 1 #index for deltamats list
  nm2 = NULL #In each pass, store the name as a comb of lag/bdy

  deltamats = list()
  for(ii in 1:nst){
    for(jj in 1:nlags){
      tm1 = transmat(tseries1,stbdyrange[[ii]],lagrange[jj])
      tm2 = transmat(tseries2,stbdyrange[[ii]],lagrange[jj])

      dt = deltatrans(tm1,tm2)
      deltamats[[ctr]] = dt$prettydelta
      nm2[ctr] = paste(names(stbdyrange)[ii],"_by_Lag",lagrange[jj],sep = "")

      ctr = ctr+1
    }
  }
  names(deltamats) = nm2

  ###Return lag range, boundary range and matrices################
  rtn = list(stbdyrange,lagrange,deltamats)
  names(rtn) = c("stbdyrange", "lagrange", "deltamats")

  return(rtn)
}
