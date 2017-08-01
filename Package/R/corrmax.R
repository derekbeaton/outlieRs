### This function should be further generalized so it can take in more compact structures, e.g., @cov and @center
  ## however, with the new version of outlieRs forthcoming, we don't retain all the copies of these types of structures, we only need to deal with X.

## For Roxygen2
#' @export
#'
#' @param population.data should be a target data set, e.g., the original raw data but only for the samples within the robust structure
#' @param sample.data should be a reference data set, e.g., the original data.

corrmax <- function(population.data,sample.data){

  return((scale(sample.data,center=colMeans(population.data),scale=apply(population.data,2,sd)) %*% sqrtMat(cor(population.data),n = -1/2))^2)

}

