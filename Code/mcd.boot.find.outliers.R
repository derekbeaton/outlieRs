mcd.boot.find.outliers <- function(X,center=F,scale=F,mcd.search.iters=1000,mcd.alpha=.75,boot.search.iters=1000,full.cloud.thresh=.99,frequency.thresh=NA){

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
