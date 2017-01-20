\name{getRMDO}
\alias{getRMDO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Mahalanobis Distance Outlyingness (RMDO) and Variable Contribution
}
\description{
Uses CovMcd() to identify the cluster with the minimum covariance determinant, and calculates the robust mean and covariance.  Then calculates the contribution of each variable to the Robust Mahalanobis Distance Outlyingness (RMDO) measure of each observation.
}
\usage{
getRMDO(dat, nsamp = "best", myalpha = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
data.frame from which the robust  mean and covariance will be identified.  Assumes complete data and variables have been centred and standardized, if necessary.
}
  \item{nsamp}{
the nsamp variable for CovMcd.
}
  \item{myalpha}{
the proportion of observations to include when searching for the minimum covariance determinant.
}
}
\details{
RMDO = 1 - (1/(1+sqrt(MD))) ### This needs to be double checked.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%% list(RMDO=rmdo, mcdCenter=x1, mcdCov=cn, initDat=dat, myalpha=myalpha, W2=w2)
  \item{RMDO}{
the RMDO value for each observation.
}
  \item{mcdCenter}{
the mean for each variable given the robust cluster of observations.
}
  \item{mcdCov}{
the covariance matrix for each variable given the robust cluster of observations.
}
  \item{initDat}{
the initial data supplied to the function.
}
  \item{myalpha}{
the alpha used in CovMcd, as supplied to the function by the user.
}
  \item{W2}{
the contribution of each variable to the observation's RMDO.  That is, an n x p data frame where higher values indicate a larger contribution.
}
%% ...
}
\references{
Garthwaite, P.H. & Koch, I. (2016).  Evaluating the contributions of individual variables to a quadratic form.  Australian & New Zealand Journal of Statistics, 58 (1): 99-119.
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
function (dat, nsamp = "best", myalpha = 0.5) 
{
    mcd1 <- CovMcd(dat, alpha = myalpha, nsamp = nsamp)
    x1 <- mcd1@center
    cn <- mcd1@cov
    md <- apply(dat, 1, function(i) {
        t(i - x1) \%*\% solve(cn) \%*\% (i - x1)
    })
    rmde <- 1 + sqrt(md)
    rmde1 <- 1/rmde
    rmdo <- 1 - rmde1
    smpsd <- sqrt(diag(cn))
    d <- diag(1/smpsd)
    dsd <- d \%*\% cn \%*\% d
    trans1 <- sweep(dat, 2, x1)
    corrmax <- sqrtMat(dsd, -1/2) \%*\% d
    contrib <- apply(trans1, 1, function(i) {
        corrmax \%*\% as.matrix(i)
    })
    w2 <- t(contrib^2)
    colnames(w2) <- colnames(dat)
    return(list(RMDO = rmdo, mcdCenter = x1, mcdCov = cn, initDat = dat, 
        myalpha = myalpha, W2 = w2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
