\name{idOutliers}
\alias{idOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying MCD-based Outliers
}
\description{
idOutliers identifies outliers from getRMDO output based on multivariate simulations of the robust mean and covariance.
}
\usage{
idOutliers(obsOV, rmdo_nsamp = 100, nsamps = 100, nreps = 100, thisquant = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsOV}{
output from getRMDO()
}
  \item{rmdo_nsamp}{
number of iterations to converge in simulated fast-MCD.  Default is 500.  "best" is optimal but would increase computation time.
}
  \item{nsamps}{
number of observations to simulate within each repetition.
}
  \item{nreps}{
number of repetitions.
}
  \item{thisquant}{
the quantile to use as the outlier cut-off.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{cutoff }{
  the cut-off based on averaging the quantiles of each repetition.
  }
  \item{cutoff_all }{
  the quantile form each repetition.
  }
%% ...
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
function (obsOV, rmdo_nsamp = 500, nsamps = 250, nreps = 100,
    thisquant = 0.95)
{
    mcdMn <- obsOV[["mcdCenter"]]
    mcdCov <- obsOV[["mcdCov"]]
    myalpha <- obsOV[["myalpha"]]
    outquants <- c()
    for (i in 1:nreps) {
        tsamp <- mvrnorm(nsamps, mcdMn, mcdCov)
        trmdo <- getRMDO(tsamp, nsamp = rmdo_nsamp, myalpha = myalpha)
        outquants[i] <- quantile(trmdo[["RMDO"]], thisquant)
    }
    simquant <- mean(outquants)
    return(list(cutoff = simquant, cutoff_all = outquants))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
