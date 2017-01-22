pca.mahal.dist <- function(X,center=T,scale=T,k=0){
	
	X <- expo.scale(X,center=center,scale=scale)
	res <- pickSVD(X,F,"svd",k=k)
	return(rowSums(res$u^2))
	
}