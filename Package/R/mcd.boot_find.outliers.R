#' Find outliers with the MCD+Bootstrap PCA approach
#' 
#' @param X Data matrix
#' @param center a boolean, vector, or string. See ?expo.scale for details.
#' @param scale a boolean, vector, or string. See ?expo.scale for details.
#' @param mcd.search.iters numeric. The number of maximum iterations to perform MCD search
#' @param mcd.alpha numeric. Threshold for proportion of MCD-based outliers during MCD search.
#' @param boot.search.iters numeric. The number of iterations to perform bootstrapping
#' @param full.cloud.thresh numeric. Threshold for proportion of values from bootstrapping to be considered outliers
#' @param frequency.thresh numeric. numeric. An inverse threshold for bootstrapping. Typically 1-full.cloud.thresh; typically unused (for now)
#' 
#' @return matrix with 8 columns. The first 4 columns are for Mahalanobis distance outlier estimates, and the last 4 columns are for Chi-squared distance outlier estimates. 
#'          the 'Dist' columns provide the distance of each item, the 'boot dist' columns provide the average of bootstrap distances
#'          the 'out %' column provides a likelihood of outlierness where the 'Outlier' columns are a simple 1 (yes) and 0 (no) of whether an item is an outlier based on the 'out %' value (and the threshold parameters input by the user)
#'          
#'
#' @examples
#' require(ExPosition)
#' data(beer.tasting.notes)
#' the.data <- expo.scale(beer.tasting.notes$data)
#' the.data_corrected_for_ABV <- apply(the.data,2, function(x){ resid(lm(x~beer.tasting.notes$sup.data[,"ABV"])) } )
#' beer_MCD.boot_full.combined_95 <- mcd.boot_find.outliers(the.data,center=F,scale=F)
#' beer.corrected_MCD.boot_full.combined_95 <- mcd.boot_find.outliers(the.data_corrected_for_ABV,center=F,scale=F)

mcd.boot_find.outliers <- function(X,center=F,scale=F,mcd.search.iters=1000,mcd.alpha=.75,boot.search.iters=1000,full.cloud.thresh=.99,frequency.thresh=NA){

	pca.res <- epPCA(X,center=center,scale=scale,graphs=F)
	mcd.boot.res <- mcd.boot.pipeline(X,center=center,scale=scale,mcd.search.iters= mcd.search.iters, mcd.alpha= mcd.alpha,boot.search.iters= boot.search.iters)
	
	
	mahal.outlier.computations <- mcd.boot.outliers(mcd.boot.res$boot.mahal, full.cloud.thresh= full.cloud.thresh,frequency.thresh= frequency.thresh,dims=dim(X),mcd.thresh=mcd.alpha)
	chi2.outlier.computations <- mcd.boot.outliers(mcd.boot.res$boot.chi2, full.cloud.thresh= full.cloud.thresh,frequency.thresh= frequency.thresh,dims=dim(X),mcd.thresh=mcd.alpha)
	
	mahal.outliers <- names(which(mahal.outlier.computations$thresh.boot.outliers))
	chi2.outliers <- names(which(chi2.outlier.computations$thresh.boot.outliers))
	
	
	outliers_full <- matrix(NA,nrow(X),8)
	rownames(outliers_full) <- rownames(X)
	colnames(outliers_full) <- c("M Dist","M mean Boot Dist","M out %","M Outlier","Chi2 Dist","Chi2 mean Boot Dist","Chi2 out %","Chi2 Outlier")
	
	
	outliers_full[names(rowSums(pca.res$ExPosition.Data$pdq$p^2)),"M Dist"] <- rowSums(pca.res$ExPosition.Data$pdq$p^2)
	outliers_full[names(rowMeans(mcd.boot.res$boot.mahal)),"M mean Boot Dist"] <- rowMeans(mcd.boot.res$boot.mahal)
	outliers_full[names(mahal.outlier.computations$freq.boot.outliers),"M out %"] <- mahal.outlier.computations$freq.boot.outliers
	outliers_full[names(mahal.outlier.computations$thresh.boot.outliers),"M Outlier"] <- mahal.outlier.computations$thresh.boot.outliers+0
	
	
	outliers_full[names(rowSums(pca.res$ExPosition.Data$fi^2)),"Chi2 Dist"] <- rowSums(pca.res$ExPosition.Data$fi^2)
	outliers_full[names(rowMeans(mcd.boot.res$boot.chi2)),"Chi2 mean Boot Dist"] <- rowMeans(mcd.boot.res$boot.chi2)
	outliers_full[names(chi2.outlier.computations$freq.boot.outliers),"Chi2 out %"] <- chi2.outlier.computations$freq.boot.outliers
	outliers_full[names(chi2.outlier.computations$thresh.boot.outliers),"Chi2 Outlier"] <- chi2.outlier.computations$thresh.boot.outliers+0
	
	
	
	return(outliers_full[order(outliers_full[,"M out %"],outliers_full[,"Chi2 out %"],decreasing=T),])
	
}
