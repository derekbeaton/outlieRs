\name{sqrtMat}
\alias{sqrtMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Square Root of a Matrix
}
\description{
Takes the square root of a matrix (if n=1/2 or some multiple).  (Not sure what happens, otherwise.  This is from stackexchange.)
}
\usage{
sqrtMat(x, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a (numeric) matrix
}
  \item{n}{
the exponent
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
function (x, n) 
{
    with(eigen(x), vectors \%*\% (values^n * t(vectors)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
