#Function to calc mean first passage time to one state from all others

mfpt = function(tm, deststate){

  ###############Check formatting#################################


  #################################################################

  #################################################################
  #Calcualte mean first passage time to reach deststate from all other states
  nn = nrow(tm)

  QQ = tm[-deststate,-deststate]
  FF = tt1 = solve(diag(nn-1)-QQ)
  mfpt = FF%*%rep(1,nn-1)

  nm = 1:nn
  nm = nm[-deststate]

  names(mfpt) = nm

  ##############Return vectors##########################
  return(mfpt)
}
