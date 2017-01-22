pca.all.dists <- function(X,center=T,scale=T,k=0,comps=2){

	if(comps <= 0){
		warning('comps must be greater than 0')
		comps <- 2
	}
			
	X <- expo.scale(X,center=center,scale=scale)
	svd.res <- pickSVD(X,F,"svd",k=k)
	
	if(comps >= length(svd.res$d)){
		warning('comps must be smaller than k: Setting comps to k-1')
		comps <- length(svd.res$d)-1
	}

		
	
	mahal.dist <- rowSums(svd.res$u^2)
	row.scores <- svd.res$u * matrix(svd.res$d,nrow(svd.res$u),ncol(svd.res$u),byrow=T)
	chi.dist <- rowSums(row.scores^2)

	dat.reduced <- (svd.res$u[,1:comps] * matrix(svd.res$d,comps,comps,byrow=T)) %*% t(svd.res$v[,1:comps])
	ortho.dist <- orthogonal.dist(X,dat.reduced)
	
	return(list(mahal.dist= mahal.dist, chi.dist= chi.dist, ortho.dist= ortho.dist))
	
}