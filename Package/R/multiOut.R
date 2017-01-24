#####################################################################################################
## Identifying Multivariate Outliers Function ##
#####################################################################################################

multiOut <- function(
  dat, # the data frame, excluding identifiers
  rmdo_alpha = 0.8, # see getRMDO -- myalpha
  rmdo_nsamp = "best", # option for getRMDO - number of iterations required to find robust parameters
  nsims = 100, # then number of total simuations (so nsims/iters per iteration)
  nsamps = 250, # the number of points to simulate in process of determining cut-off for outlier
  alpha=0.01 # error to allow in CIs for quantile cut-offs
){
  
  # Giving column and row names if none are previously present
  if(is.null(colnames(dat))){
    colnames(dat) <- paste("Col",1:ncol(dat),sep="")
  } 
  if(is.null(rownames(dat))){
    rownames(dat) <- paste("Row",1:nrow(dat),sep="")
  } 
  
  # Calling getRMDO to calculate outlyingness
  rmdos <- getRMDO(dat,nsamp=rmdo_nsamp,myalpha=rmdo_alpha)
  
  # Output of getRMDO
  rmdoCalc <- rmdos[["RMDO"]]
  centers <- rmdos[["mcdCenter"]]
  scatters <- data.frame(rmdos[["mcdCov"]])
  contVar <- data.frame(rmdos[["W2"]])
  contVar$ID <- rownames(contVar)
  
  # Simulating to identify "normal range"
  idOut <- idOutliers(obsOV=rmdos,nsamps=nsamps,nreps=nsims,thisquant = 1-alpha)
  outCut <- idOut[["cutoff"]]
  
  # Reformatting contVar to long form for ease of interpretation
  contVar1 <- suppressMessages(melt(contVar))
  
  # Making final lists of outliers
  outs <- data.frame(ID=rownames(dat),RMDO=rmdoCalc,Outlier=ifelse(rmdoCalc >= outCut,1,0))
  outs1 <- outs[order(outs$RMDO,decreasing = T),]
  
  basicOut <- outs[,c("ID","Outlier")]
  colnames(basicOut)[2] <- "mcd_outlier"
  
  return(list(outlier_decision=basicOut,
              hidden_detail=list(initDat=as.data.frame(dat),RMDO=outs1, RMDOalpha=rmdo_alpha, Contribs=contVar1, CutOff=outCut)))
} # end of multiOut function
