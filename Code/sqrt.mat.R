##a function specifically for sqrt'ing a matrix
##stolen from http://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-


###this could be sped up at the bottom.
sqrt.mat <- function(X,tol=sqrt(.Machine$double.eps)){
	
	isDiagonal.matrix <- function(X){
		if(is.null(dim(X))){
			stop("sqrt.mat: X is not a matrix.")
		}
		return(all(X[lower.tri(X)] == 0, X[upper.tri(X)] == 0))
		
	}
	
	##first, test if it is symmetric in size; and other symmetric properties.
	if(!isSymmetric.matrix(X,tol=tol)){
		warning("sqrt.mat: Weight/Mass matrix is not symmetric; using SVD for sqrt of rectangular matrix.")
		return(svd.sqrt.mat(X))
		###I can just do the SVD one here...
	}
			
	##test if it is a diagonal matrix -- then just sqrt it and send back.
	if(isDiagonal.matrix(X)){
		return(  diag(sqrt(diag(X))) )
	}else{
		
		A <- eigen(X)
		##test for complex values:
		if(any(unlist(lapply(A$values,is.complex)))){
			stop("sqrt.mat: Weight/Mass matrix not positive definitive. Eigenvalues are complex.")
		}
		
		##change values below tolerance		
		A$values[which(A$values < tol)] <- 0

		##first, test if positive definite
		if( sum(A$values < 0 )>0 ){
			stop("sqrt.mat: Weight/Mass matrix not positive definite. Some eigenvalues are less than 0")	
		}else{		
			return((A$vectors * matrix(sqrt(A$values),nrow(A$vectors),ncol(A$vectors),byrow=T)) %*% t(A$vectors))
		}
	}
}
