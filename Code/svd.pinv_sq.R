###Am using this in place of ginv() for square symmetric matrices. ginv() returns some rounding errors and makes the eigenvalues imaginary. 
svd.pinv_sq <- function(X,tol=sqrt(.Machine$double.eps)){
	
	###I am requiring this to be square.
	if( !isSymmetric.matrix(X, tol= tol) ){
		stop("svd.pinv_sq: Matrix not symmetric")
	}
	res <- svd(X)
	drop.vecs <- which(res$d < tol)
	if(length(drop.vecs) > 0){
		res$d <- res$d[-c(drop.vecs)]
		res$u <- res$u[,-c(drop.vecs)]
			res$u[abs(res$u) < tol] <- 0
	}
	return( (res$u * matrix(1/res$d,nrow(res$u),ncol(res$u),byrow=T)) %*% t(res$u) )
}