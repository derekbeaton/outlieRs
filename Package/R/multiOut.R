#' @export
#'
#' @description Calls getRMDO() and idOutliers() to calculate the robust parameters, RMDO, and identify the cutoff for outliers.  Then, compares the RMDO for each observation with the outlier cutoff.  It also reformats the contribution matrix.
#' @title multivariate outliers (via MCD)
#'
#' @param dat the data frame from which to identify outliers. Observation identifiers should be in the rownames. Assumes complete data and variables have been centred and standardized, if necessary.
#' @param rmdo_alpha the alpha to be used in CovMcd – the proportion of observations to keep in the robust cluster.
#' @param rmdo_samp the nsamp to be used in CovMcd – the maximum number of iterations to use for convergence of the robust cluster. "best" is recommended as "exhaustive enumeration is done, as long as the number of trials does not exceed 5000." See also ?\code{CovMcd}
#' @param nsims number of repetitions of simulating multivariate normal in idOutliers.
#' @param nsamps number of observations to simulate within each repetition for idOutliers.
#' @param alpha the error to allow for quantile cut-offs. 1-alpha -> thisquant in idOutliers.
#'
#' @return
#' A list with two elements: outlier_decision and hidden_detail.
#' \item{outlier_decision}{ a data.frame indicating the observations that are outliers: 1 = outlier; 0 = not an outlier.}
#' hidden_detail is another list with the additional information:
#' \item{initDat}{the initial data supplied to the function.}
#' \item{RMDO}{data.frame of the outlier results in decreasing order of most-least outlying.Outlier=1 implies the observation is an outlier, while Outlier=0 implies otherwise.}
#' \item{RMDOalpha}{the alpha used in CovMcd, as supplied to the function by the user.}
#' \item{Contribs}{data.frame of the contribution for each observation-variable pair.}
#' \item{Cutoff}{the cut-off based on averaging the quantiles of each repetition.}
#'
#' @examples
#' data(beer.tasting.notes)
#' the.data <- expo.scale(beer.tasting.notes$data)
#' the.data_corrected_for_ABV <- apply(the.data,2, function(x){ resid(lm(x~beer.tasting.notes$sup.data[,"ABV"])) } )
#' beer_mcd <- multiOut(the.data)
#' beer.corrected_mcd <- multiOut(the.data_corrected_for_ABV)

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
              hidden_detail=list(initDat=as.data.frame(dat),RMDO=outs1, RMDOalpha=rmdo_alpha, Contribs=contVar1, corrMat=rmdos[["corrMat"]], CutOff=outCut)))
} # end of multiOut function
