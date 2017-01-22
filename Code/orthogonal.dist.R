orthogonal.dist <- function(X,Y){
	
	if( nrow(X) == nrow(Y) & ncol(X) == ncol(Y) ){
		return( rowSums((X - Y)^2) )
	}else{
		stop("X and Y do not have the same dimensions")
	}
	
}