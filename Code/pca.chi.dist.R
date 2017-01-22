pca.mahal.dist <- function(X,center=T,scale=T,k=0){

	X <- expo.scale(X,center=center,scale=scale)
	res <- pickSVD(X,k=k)
	return( rowSums( (res$u * matrix(res$d,nrow(res$u),length(res$d),byrow=T))^2 ) )
	
}