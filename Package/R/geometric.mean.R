### This function should be further generalized so it can take in more compact structures, e.g., @cov and @center
## however, with the new version of outlieRs forthcoming, we don't retain all the copies of these types of structures, we only need to deal with X.

## For Roxygen2
#' @export
#'
#' @title geometric.mean -- STOLEN FROM PSYCH PACKAGE.
#'
#' @param x vector
#' @param na.rm flag to ignore NAs
#' 
#' @author William Revelle
#' 
#' @seealso \code{\link{psych}}
#' 

geometric.mean <- function (x, na.rm = TRUE) 
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}