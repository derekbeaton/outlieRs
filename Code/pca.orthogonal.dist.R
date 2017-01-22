pca.orthogonal.dist <- function(X,center=T,scale=T,k=5){
	
	
	X <- expo.scale(X,center=center,scale=scale)
	svd.res <- pickSVD(X,F,"svd",k=k)
	dat.reduced <- (svd.res$u * matrix(svd.res$d,nrow(svd.res$u),ncol(svd.res$u),byrow=T)) %*% t(svd.res$v)
	return( orthogonal.dist(X,dat.reduced) )
	
}