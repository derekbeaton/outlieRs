#####################################################################################################
## Calculating the Robust MCD-based Mahalanobis Distance Outlyingness measure ##
#####################################################################################################
getRMDO <- function(
  dat,
  nsamp="best", # the nsamp variable for CovMcd
  myalpha = 0.5 # alpha for proportion of obs to include in MCD-MD
){
  # Using fast-MCD to calculate robust mean and covariance
  mcd1 <- CovMcd(dat,alpha=myalpha,nsamp=nsamp) 
  x1 <- mcd1@center
  cn <- mcd1@cov
  
  # Calculating the Mahalanobis distance (MD), and then robust MD outlyingness (RMDO) based on robust mean and covariance
  md <- apply(dat,1,function(i){t(i-x1)%*%solve(cn)%*%(i-x1)})
  rmde <- 1 + sqrt(md)
  rmde1 <- 1/rmde
  rmdo <- 1-rmde1 # outlyingness values
  #   
  #   used <- rownames(dat)[mcd1@best] # the subjects who were included in the MCD calculations
  #   
  # Based on 2016 paper: Contributions to Quadratic Form
  smpsd <- sqrt(diag(cn))
  d <- diag(1/smpsd)
  dsd <- d %*% cn %*% d
  
  trans1 <- sweep(dat,2,x1)
  corrmax <- sqrtMat(dsd,-1/2) %*% d
  contrib <- apply(trans1,1,function(i){corrmax %*% as.matrix(i)})
  w2 <- t(contrib**2)
  colnames(w2) <- colnames(dat)
  
  return(list(RMDO=rmdo, mcdCenter=x1, mcdCov=cn, initDat=dat, myalpha=myalpha, W2=w2))
  # mcdUsed=used, DSD=dsd, CorrMax=corrmax, 
}

