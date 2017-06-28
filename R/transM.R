transM = function(x,binbounds,lag){##Function to find transition prob matrix
  ##Create lagged variable
  nbds = length(binbounds); nbins = nbds-1
  
  tmp = x
  tmp$X = 1:nrow(tmp)
  tmp$lag = ifelse(tmp$X+lag<=nrow(tmp), tmp$MassAve[tmp$X+lag], NA)
  
  tmp$sourcebin = findInterval(tmp$MassAve,binbounds)
  tmp$destbin = findInterval(tmp$lag,binbounds)
  
  tM =  matrix(0, nrow = nbins, ncol = nbins)
  for(ii in 1:nbins){
    for(jj in 1:nbins){
      tM[ii,jj] = length(which(tmp$sourcebin == ii & tmp$destbin == jj))
    }
  }
  tM.rs = rowSums(tM)
  tM.prob = tM/tM.rs
  
  return(list(tM,tM.prob))
}