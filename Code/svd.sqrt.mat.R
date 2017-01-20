##an alternative square root function for rectangular matrices (if needed).
svd.sqrt.mat <- function(X,tol=sqrt(.Machine$double.eps)){
	# if( !isSymmetric.matrix(X) ){
		# stop("svd.sqrt.mat: Matrix not symmetric")
	# }	
	res <- svd(X)
	drop.vecs <- which(res$d < tol )
	if(length(drop.vecs) > 0){
		res$d <- res$d[-c(drop.vecs)]
		res$u <- res$u[,-c(drop.vecs)]
			res$u[abs(res$u) < tol] <- 0
		res$v <- res$v[,-c(drop.vecs)]		
			res$v[abs(res$v) < tol] <- 0
	}
	return( (res$u * matrix(sqrt(res$d),nrow(res$u),ncol(res$u),byrow=T)) %*% t(res$v) )
}