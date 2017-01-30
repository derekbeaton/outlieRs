# Find outliers with the Candes Robust PCA approach
#' @export

### DO THESE LATER.
# @description Uses CovMcd() to identify the cluster with the minimum covariance determinant, and calculates the robust mean and covariance.  Then calculates the contribution of each variable to the Robust Mahalanobis Distance Outlyingness (RMDO) measure of each observation.
#
# @param dat data.frame from which the robust  mean and covariance will be identified.  Assumes complete data and variables have been centred and standardized, if necessary.
# @param nsamp the nsamp variable for CovMcd.
# @param myalpha the proportion of observations to include when searching for the minimum covariance determinant.
#
# @details RMDO = 1 - (1/(1+sqrt(MD))) ### This needs to be double checked.
#

getRMDO <- function(dat,nsamp="best", myalpha = 0.5 ){

	### DB made this change!
	#### Here I am moving the sqrtMat() function from its own file to a "private" function within getRMDO.
	#### the reason is that it is a function used here, that we don't necessarily want the users to access.
	#####################################################################################################

	## a private (to getRMDO) function: sqrtMat
	## function to calculate cov^(-1/2) ##
	#####################################################################################################
	sqrtMat <- function(x, n){
	  with(eigen(x), vectors %*% (values^n * t(vectors)))
	}


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
  d <- diag(1/smpsd)      ## DB NOTE: a fraction can make a diag() panic.
  dsd <- d %*% cn %*% d   ## DB Q: Is this converting cov to cor?

  trans1 <- sweep(dat,2,x1) ## DB Q: Center data by robust mean?
  corrmax <- sqrtMat(dsd,-1/2) %*% d ## DB Q: inverse sqrt of cor times data?
  colnames(corrmax) <- rownames(corrmax) <- colnames(dat)
  contrib <- apply(trans1,1,function(i){corrmax %*% as.matrix(i)}) # DB Q: data times correlation?
  w2 <- t(contrib**2) ## DB Q: squared "contributions"?
  colnames(w2) <- colnames(dat)

  return(list(RMDO=rmdo, mcdCenter=x1, mcdCov=cn, initDat=dat, myalpha=myalpha, W2=w2,corrMat=corrmax))
  # mcdUsed=used, DSD=dsd, CorrMax=corrmax,
}

