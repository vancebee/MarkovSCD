#Function to calculate and determine the significance between two transition matrices

deltatrans = function(tm1, tm2, nruns = 10000, CI = 0.95, mineffect = 0.05){

  ###############Check formatting#################################
  if(CI <=0 | CI>= 1) stop("Confidence interval must be between 0 and 1")

  ############Difference between transition matrices###############
  delta = tm2$prob-tm1$prob
  delta = ifelse(is.na(delta),0,delta)
  nn = nrow(delta)

  #################Upper/Lower Index###############################
  l = round(((1-.95)/2)*nruns); u = round((1-(1-.95)/2)*nruns) #Upper and lower index

  #################################################################
  #################################################################
  ###Inline function to generate multinomial sample from a row#####
  multibs = function(x){
    nn = length(x) #nbins+1
    if(!x[nn] == 0){
      bs.count = rmultinom(nruns,x[nn],x[1:(nn-1)]) #Generate random multinomial
      bs.prob = apply(bs.count,2,function(y) y/x[nn]) #Convert to probability
    }else matrix(0,nrow = 6, ncol = nruns) #If no entries, return LCI and UCI as 0,0
  }
  #################################################################
  #################################################################

  ###Append total # obs. in each state to prob transition matrix###
  rs1 = rowSums(tm1$raw) #total entries in each row
  boot1 = cbind(tm1$prob,rs1)
  rs2 = rowSums(tm2$raw) #total entries in each row
  boot2 = cbind(tm2$prob,rs2)

  #################################################################
  #Create list of generated probabilities. One list element per
  #state.  Each element is a matrix with column representing a
  #bootstrap run and each row having the probability for a state
  CIlist1 = plyr::alply(boot1,1,multibs)
  CIlist2 = plyr::alply(boot2,1,multibs)

  #################################################################
  #Calc differences between the proabilities on a
  #per bootstrap run basis
  deltaboot = list()
  for(ii in 1:length( CIlist1)){
    deltaboot[[ii]] = CIlist1[[ii]]-CIlist2[[ii]]
  }

  #################################################################
  #Sort all of the bootstrap differences and select CI boundaries
  #based on CI specification.  Output is a list.  Each element
  #corresponds to a row in transition matrix and consists of a list
  #of vectors, lower then upper
  srtdelta = lapply(deltaboot, function(z) t(apply(z,1,sort)))
  lu = lapply(srtdelta, function(z) list(z[,l],z[,u]))

  ################Determine if CI interval contains 0##############
  sigmat = matrix(NA,nn,nn)
  for(ii in 1:nn){
    for(jj in 1:nn){
      su = sign(lu[[ii]][[1]][jj])
      sl = sign(lu[[ii]][[2]][jj])
      sigmat[ii,jj] = ifelse(su==sl & !(su==0 & sl==0),1,0)
    }
  }

  ######Only allow significance if effect is large enough##########
  sigmat[which(abs(delta)<0.05)] = 0
  rownames(sigmat) = rownames(delta)
  colnames(sigmat) = colnames(delta)

  #################################################################
  #Create delta matrix that only prints sig values & rounds values#
  prettydelta = matrix(NA,nn,nn)
  for(ii in 1:nrow(delta)){
    for(jj in 1:ncol(delta)){
      prettydelta[ii,jj] = ifelse(sigmat[ii,jj]==1,round(delta[ii,jj],2),0)
    }
  }
  rownames(prettydelta) = rownames(delta)
  colnames(prettydelta) = colnames(delta)

  ##############Return list of matrices##########################
  rtn = list(delta,sigmat,prettydelta)
  names(rtn) = c("delta", "sigmat", "prettydelta")
  return(rtn)
}
