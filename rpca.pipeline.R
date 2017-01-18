rpca.pipeline <- function(X,S.tolerance=.Machine$double.eps * 100,
					default.lambda=1/sqrt(max(dim(X))),
					lambda.vec=(1/sqrt(min(dim((X))) / seq(1,10,1))),
					sparsity.threshold=0.05, corrupted.tolerance=0){
	
	X.dim <- dim(X)
	X.cells <- X.dim[1] * X.dim[2]
	
	print("Begin 'default' robust PCA")
	rob.res_def <- rpca(X)
	rownames(rob.res_def$S) <- rownames(rob.res_def$L) <- rownames(X)
	colnames(rob.res_def$S) <- colnames(rob.res_def$L) <- colnames(X)
	rob.res_def$S[abs(rob.res_def$S) < S.tolerance] <- 0
	print("End 'default' robust PCA")


	print("Begin 'sparse' robust PCA search")
	print(paste0("Length of search vector: ",length(lambda.vec)))
	for(i in 1:length(lambda.vec)){
		rob.res_sparse <- rpca(X,lambda=lambda.vec[i])
		corrupted.sum <- sum(abs(rob.res_sparse$S) > S.tolerance)
		print(i)
		if(	(corrupted.sum / X.cells) < (sparsity.threshold + corrupted.tolerance)	){
			break
		}
	}
	rob.res_sparse$S[abs(rob.res_sparse$S) < S.tolerance] <- 0
	rownames(rob.res_sparse$S) <- rownames(rob.res_sparse$L) <- rownames(X)
	colnames(rob.res_sparse$S) <- colnames(rob.res_sparse$L) <- colnames(X)
	print("End 'sparse' robust PCA search")	

	return(list(rpca.default= rob.res_def,rpca.sparse= rob.res_sparse))
	
}