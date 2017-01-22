##for now, it's just the SVD, but could easily be any actual SVD technique
fast.pc.dists <- function(orig.data,o.center=F,o.scale=F,target.data,t.center=F,t.scale=F,t.comps=0){
	
	if( t.comps > max(dim(target.data))){
		warning("t.comps is too high, changing to min(dim(target.data))")
		t.comps <- min(dim(target.data))
	}
	
	##theoretically, we could also use a subspace of the original, too, but I'd rather not do that within this code.
	
	orig.data <- expo.scale(orig.data,center=o.center,scale=o.scale)
	target.data <- expo.scale(target.data,center=t.center,scale=t.scale)	
	
	##for fastness, just SVD. But this has to be generalized to whatever method...
	target.svd <- pickSVD(target.data,k=t.comps)
	row.scores <- target.svd$u * matrix(target.svd$d,nrow(target.svd$u),length(target.svd$d),byrow=T)
	target.mahal <- rowSums(target.svd$u^2)
	target.chi2 <- rowSums( (row.scores)^2 )
	
	
	if( t.comps==0 | t.comps <= min(dim(target.data)) ){
		orth.dist <- rowSums((orig.data - target.data)^2)
	}else{
		orth.dist <- rowSums( (orig.data - (row.scores %*% t(target.svd$v)) )^2 )
	}
	
	return(list(target.mahal=target.mahal,target.chi2=target.chi2, OtoT_orth.dist=orth.dist))
	
	##NOTE with MCA, it will be either a hamming-like distance of the 0s/1s or it should be done on the matrix to be decomposed.
	
	
}