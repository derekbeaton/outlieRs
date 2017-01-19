rpca.outliers <- function(X,quantile.thresh=.9,S.tolerance=.Machine$double.eps * 100,
					default.lambda=1/sqrt(max(dim(X))),
					lambda.vec=(1/sqrt(min(dim((X))) / seq(1,10,1))),
					sparsity.threshold=0.05, corrupted.tolerance=0){


	rpca.pipe.res <- rpca.pipeline(X, S.tolerance= S.tolerance, default.lambda= default.lambda, lambda.vec= lambda.vec, sparsity.threshold= sparsity.threshold, corrupted.tolerance= corrupted.tolerance)
	
	S.thresh <- rpca.pipe.res$rpca.sparse$S
	S.thresh[S.thresh < 0] <- -1
	S.thresh[S.thresh > 0] <- 1
	corrupted.rows <- which(rowSums(abs(S.thresh)) > 0)	
	corrupted.cols <- which(colSums(abs(S.thresh[corrupted.rows,]))>0)
	corrupted.threshold <- S.thresh[corrupted.rows,corrupted.cols]
	
	###distances from SVDs
	X.from.S_def <- fast.pc.dists(X,target.data=rpca.pipe.res$rpca.default$S)
	X.from.L_def <- fast.pc.dists(X,target.data=rpca.pipe.res$rpca.default$L)
	L_def.from.S_def  <- fast.pc.dists(rpca.pipe.res$rpca.default$L,target.data=rpca.pipe.res$rpca.default$S)
	
	X.from.S_sparse <- fast.pc.dists(X,target.data=rpca.pipe.res$rpca.sparse$S)
	X.from.L_sparse <- fast.pc.dists(X,target.data=rpca.pipe.res$rpca.sparse$L)
	L_sparse.from.S_sparse  <- fast.pc.dists(rpca.pipe.res$rpca.sparse$L,target.data=rpca.pipe.res$rpca.sparse$S)
	
	low.rank.proportional.dist <- X.from.S_def$OtoT_orth.dist / L_def.from.S_def$OtoT_orth.dist		##anything here below 1
	low.rank.proportional.dist[is.nan(low.rank.proportional.dist)] <- 0	
	low.rank.proportional.dist <- round(low.rank.proportional.dist,digits=8)
	
		## a sufficient threshold
	sparse.low.rank.proportional.dist <- X.from.S_sparse$OtoT_orth.dist / L_sparse.from.S_sparse$OtoT_orth.dist	##anything here *not 1*
	sparse.low.rank.proportional.dist[is.nan(sparse.low.rank.proportional.dist)] <- 0	
	sparse.low.rank.proportional.dist <- round(sparse.low.rank.proportional.dist,digits=8)
			
	low.rank.outliers <- which((1/low.rank.proportional.dist) > quantile(1/low.rank.proportional.dist,probs= quantile.thresh))
	sparse.low.rank.outliers <- which((sparse.low.rank.proportional.dist) > quantile(sparse.low.rank.proportional.dist,probs= quantile.thresh))
	
	
	###PERHAPS I NEED SOMETHING BETTER THAN QUANTILES?!
	
	
	outliers <- matrix(0,nrow(X),3)
	rownames(outliers) <- rownames(X)
	colnames(outliers) <- c("Low Rank > 90% Quant","Sparse Low Rank > 90% Quant","RPCA Corrupted")
	outliers[names(low.rank.outliers),"Low Rank > 90% Quant"] <- 1
	outliers[names(sparse.low.rank.outliers),"Sparse Low Rank > 90% Quant"] <- 1
	outliers[names(corrupted.rows),"RPCA Corrupted"] <- 1
	
	
	return(list(outliers=outliers,corrupted.threshold= corrupted.threshold,
		dists=list(
				X.from.S_def_OtoT_orth.dist=X.from.S_def$OtoT_orth.dist,
				L_def.from.S_def_OtoT_orth.dist=L_def.from.S_def$OtoT_orth.dist,
				X.from.S_sparse_OtoT_orth.dist=X.from.S_sparse$OtoT_orth.dist,
				L_sparse.from.S_sparse_OtoT_orth.dist=L_sparse.from.S_sparse$OtoT_orth.dist
				)
		))

	
}