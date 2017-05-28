# library(InPosition)
# library(rrcov)
# library(doParallel)

mcd.boot.pipeline <- function(X, center=T,scale=T, mcd.search.iters=2500, mcd.alpha=.75, boot.search.iters=1000){
	
	print(paste0("MCD Iters: ", mcd.search.iters))
	print(paste0("Boot Iters: ", boot.search.iters))
	print(mcd.alpha)

	
	print('Start initial PCA')	
	fixed.res <- epPCA(X,center=center,scale=scale,graphs=F)
	no.varX <- fixed.res$ExPosition.Data$pdq$p %*% t(fixed.res$ExPosition.Data$pdq$q)
	print('End initial PCA')			
		
	##just the MCD-like approach for now.
	sub.samp.size <- min(h.alpha.n(mcd.alpha,nrow(X),ncol(X)),nrow(X)-1)
	sub.samps <- foreach(i=1:mcd.search.iters,.combine=rbind) %dopar% sort(sample(nrow(X), sub.samp.size))
	sub.samps <- unique(sub.samps)
	sub.iters <- nrow(sub.samps)
	min.log.det <- sum(log(fixed.res$ExPosition.Data$eigs/ (nrow(X)-1)))
	min.mcd.samps <- c()
	print('Start MCD search')	## I can make this better. I need to focus on the WIRE-CS paper.
								## but I think we can consider a pseudo-split-half approach later.
	pb <- txtProgressBar(1, sub.iters, 1, style = 1)
	for(i in 1:sub.iters){
		sub.X <- X[sub.samps[i,],]
		sub.eigs <- epPCA(sub.X,center=center,scale=scale,graphs=F)$ExPosition.Data$eigs
		this.log.det <- sum(log(sub.eigs/ (nrow(sub.X)-1)))
		if(this.log.det < min.log.det){
			min.log.det <- this.log.det
			min.mcd.samps <- sub.samps[i,]
		}
		setTxtProgressBar(pb, i)
	}
	close(pb)
	print('End MCD search')	
	sub.pca.mcd <- epPCA(X[min.mcd.samps,],center=center,scale=scale,graphs=F)
	
	boot.samps_mcd <- foreach(i=1: boot.search.iters,.combine=rbind) %dopar% sort(sample(1:length(min.mcd.samps),replace=T))
	boot.samps_mcd <- unique(boot.samps_mcd)
	boot.iters <- nrow(boot.samps_mcd)
	boot.FJs_mcd <- array(NA, dim=c( ncol(X), length(sub.pca.mcd$ExPosition.Data$eigs), boot.iters) )
	boot.FIs_mcd <- array(NA, dim=c( nrow(X), length(sub.pca.mcd$ExPosition.Data$eigs), boot.iters) )
	rownames(boot.FJs_mcd) <- colnames(X)
	rownames(boot.FIs_mcd) <- rownames(X)	
	
	print('Start Bootstrap')
	pb <- txtProgressBar(1, boot.iters, 1, style = 1)
	for(i in 1:boot.iters){
		boot.X_mcd <- X[min.mcd.samps[boot.samps_mcd[i,]],]
		boot.ps_mcd <- sub.pca.mcd$ExPosition.Data$pdq$p[boot.samps_mcd[i,],] 
		boot.FJs_mcd[,,i] <- t(expo.scale(boot.X_mcd,center = sub.pca.mcd$ExPosition.Data$center, scale = sub.pca.mcd$ExPosition.Data$scale)) %*% boot.ps_mcd
		boot.FIs_mcd[,,i] <- no.varX %*% boot.FJs_mcd[,,i]
		setTxtProgressBar(pb, i)
	}
	close(pb)
	print('End Bootstrap')
	
	rownames(boot.FIs_mcd) <- rownames(X)
	rownames(boot.FJs_mcd) <- colnames(X)	
	
	print('Start: Compute Outliers')	
	boot.mahal.mcd <- apply(boot.FIs_mcd,3,function(x) {rowSums(expo.scale(x,center=T,scale="SS1")^2)})
	boot.chi2.mcd <- apply(boot.FIs_mcd,3,function(x) {rowSums(x^2)})
	print('End: Compute Outliers')
	

	return(list(	mcd.list = min.mcd.samps,
					mcd.log.det = min.log.det,
					boot.mahal = boot.mahal.mcd,
					boot.chi2 = boot.chi2.mcd))
	
}
