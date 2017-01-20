###not using this as of now, but can later.
svd.pinv <- function(X,tol=sqrt(.Machine$double.eps)){
	res <- svd(X)
	drop.vecs <- which(res$d < tol )
	if(length(drop.vecs) > 0){
		res$d <- res$d[-c(drop.vecs)]
		res$u <- res$u[,-c(drop.vecs)]
			res$u[abs(res$u) < tol] <- 0
		res$v <- res$v[,-c(drop.vecs)]		
			res$v[abs(res$v) < tol] <- 0
	}	
	return( (res$v * matrix(1/res$d,nrow(res$v),ncol(res$v),byrow=T)) %*% t(res$u) )
}
