\name{multiOut}
\alias{multiOut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Outliers with MCD
}
\description{
Calls getRMDO() and idOutliers() to calculate the robust parameters, RMDO, and identify the cutoff for outliers.  Then, compares the RMDO for each observation with the outlier cutoff.  It also reformats the contribution matrix.
}
\usage{
multiOut(dat, rmdo_alpha = 0.8, rmdo_nsamp = "best", nsims = 100, nsamps = 250, alpha = 0.01)
}


\arguments{
  \item{dat}{
the data frame from which to identify outliers.  Observation identifiers should be in the rownames.  Assumes complete data and variables have been centred and standardized, if necessary.
}
  \item{rmdo_alpha}{
the alpha to be used in CovMcd -- the proportion of observations to keep in the robust cluster.
}
  \item{rmdo_nsamp}{
the nsamp to be used in CovMcd -- the maximum number of iterations to use for convergence of the robust cluster.  "best" is recommended as "exhaustive enumeration is done, as long as the number of trials does not exceed 5000."
}
  \item{nsims}{
number of repetitions of simulating multivariate normal in idOutliers.
}
  \item{nsamps}{
number of observations to simulate within each repetition for idOutliers.
}
  \item{alpha}{
the error to allow for quantile cut-offs. 1-alpha -> thisquant in idOutliers.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%% return(list(initDat=as.data.frame(dat),RMDO=outs1, RMDOalpha=rmdo_alpha, Contribs=contVar1, CutOff=outCut))
A list with two elements: outlier_decision and hidden_detail.
  \item{outlier_decision}{
a data.frame indicating the observations that are outliers: 1 = outlier; 0 = not an outlier.
}

hidden_detail is another list with the additional information:

  \item{initDat}{
the initial data supplied to the function.
}
  \item{RMDO}{
data.frame of the outlier results in decreasing order of most-least outlying.Outlier=1 implies the observation is an outlier, while Outlier=0 implies otherwise.
}
  \item{RMDOalpha}{
the alpha used in CovMcd, as supplied to the function by the user.
}
  \item{Contribs}{
data.frame of the contribution for each observation-variable pair.
}
  \item{Cutoff}{
  the cut-off based on averaging the quantiles of each repetition.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, rmdo_alpha = 0.8, rmdo_nsamp = "best", nsims = 100,
    nsamps = 250, alpha = 0.01)
{
    if (is.null(colnames(dat))) {
        colnames(dat) <- paste("Col", 1:ncol(dat), sep = "")
    }
    if (is.null(rownames(dat))) {
        rownames(dat) <- paste("Row", 1:nrow(dat), sep = "")
    }
    rmdos <- getRMDO(dat, nsamp = rmdo_nsamp, myalpha = rmdo_alpha)
    rmdoCalc <- rmdos[["RMDO"]]
    centers <- rmdos[["mcdCenter"]]
    scatters <- data.frame(rmdos[["mcdCov"]])
    contVar <- data.frame(rmdos[["W2"]])
    contVar$ID <- rownames(contVar)
    idOut <- idOutliers(obsOV = rmdos, nsamps = nsamps, nreps = nsims,
        thisquant = 1 - alpha)
    outCut <- idOut[["cutoff"]]
    contVar1 <- suppressMessages(melt(contVar))
    outs <- data.frame(ID = rownames(dat), RMDO = rmdoCalc, Outlier = ifelse(rmdoCalc >=
        outCut, 1, 0))
    outs1 <- outs[order(outs$RMDO, decreasing = T), ]
    return(list(initDat = as.data.frame(dat), RMDO = outs1, RMDOalpha = rmdo_alpha,
        Contribs = contVar1, CutOff = outCut))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
