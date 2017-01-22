#' Find outliers with the Candes Robust PCA approach
#' 
#' @param X Data matrix
#' @param center a boolean, vector, or string. See ?expo.scale for details.
#' @param scale a boolean, vector, or string. See ?expo.scale for details.
#' @param quantile.thresh numeric. 
#' @param S.tolerance numeric. A tolerance value to threshold extremely small numbers in the S(parse) matrix returned by rpca()
#' @param default.lambda numeric. The default lambda value for rpca, which is used for the thresholding procedure
#' @param lambda.vec vector (of numerics). A sequence of lambda values used during the sparsity.threshold search; if sparsity.threshold is not met, rpca sparsity search will stop after it goes through all of these values.
#' @param sparsity.threshold numeric. A threshold of how sparse the S(parse) matrix should be. lambda.vec is used to perform this search and when sparsity is satisified from that vector, rpca returns its results. 
#' @param corrupted.tolerance numeric. A tolerance factor to adjust for \code{sparsity.threshold} to allow for 'close enough' +/- tolerance.
#' 
#' @return list with many 4 items: 
#'          $outliers -- a matrix with 4 columns. 
#'              The first column is 1 (yes) or 0 (no) based on the quantile (based on \code{quantile.thresh}) cutoff based on the proportion of the orthogonal distances: \code{orthogonal.dist(X,S) / orthogonal.dist(L,S)} based on the 'default' RPCA (which uses \code{default.lambda}) 
#'              The second column is 1 (yes) or 0 (no) based on the quantile (based on \code{quantile.thresh}) cutoff based on the proportion of the orthogonal distances: \code{orthogonal.dist(X,S) / orthogonal.dist(L,S)} based on the S(parse) and L(ow rank) matrices discovered in the search process (via \code{lambda.vec} and \code{sparsity.threshold})
#'              The third column is 1 (yes) or 0 (no) based on whether the default RPCA found any data corruptions for that row
#'              The fourth column is 1 (yes) or 0 (no) based on whether the sparsity search RPCA found any data corruptions for that row
#' 
#'          $def.corrupted.threshold (from the default RPCA) a thresholded matrix of only the corrupted rows that also indicates which columns were corrupted
#'          
#'          $corrupted.threshold (from the sparsity search) a thresholded matrix of only the corrupted rows that also indicates which columns were corrupted
#'          
#'          $dists - a list of the orthogonal distances used to compute outlier estimates in $outliers:
#'              X.from.S_def_OtoT_orth.dist -- \code{orthogonal.distance(X,S)} from default RPCA
#'              L_def.from.S_def_OtoT_orth.dist -- \code{orthogonal.distance(L,S)} from default RPCA              
#'              X.from.S_sparse_OtoT_orth.dist -- \code{orthogonal.distance(X,S)} from sparsity search RPCA
#'              L_sparse.from.S_sparse_OtoT_orth.dist -- \code{orthogonal.distance(L,S)} from sparsity search RPCA                            
#'              
#' @examples 
#' require(ExPosition)
#' data(beer.tasting.notes)
#' the.data <- expo.scale(beer.tasting.notes$data)
#' the.data_corrected_for_ABV <- apply(the.data,2, function(x){ resid(lm(x~beer.tasting.notes$sup.data[,"ABV"])) } )
#' beer_rpca <- rpca_find.outliers(the.data)
#' beer.corrected_rpca <- rpca_find.outliers(the.data_corrected_for_ABV)



rpca_find.outliers <- function(X, center=F, scale= F,quantile.thresh=.9,S.tolerance=.Machine$double.eps * 100,
					default.lambda=1/sqrt(max(dim(X))),
					lambda.vec=(1/sqrt(min(dim((X))) / seq(1,10,1))),
					sparsity.threshold=0.05, corrupted.tolerance=0){


  if( length(lambda.vec) < 1 | is.null(lambda.vec)){
    lambda.vec <- default.lambda
  }
  
	rpca.pipe.res <- rpca.pipeline(X, center=center, scale=scale, S.tolerance= S.tolerance, default.lambda= default.lambda, lambda.vec= lambda.vec, sparsity.threshold= sparsity.threshold, corrupted.tolerance= corrupted.tolerance)
	
	S.thresh <- rpca.pipe.res$rpca.sparse$S
	S.thresh[S.thresh < 0] <- -1
	S.thresh[S.thresh > 0] <- 1
	corrupted.rows <- which(rowSums(abs(S.thresh)) > 0)	
	corrupted.cols <- which(colSums(abs(S.thresh[corrupted.rows,]))>0)
	corrupted.threshold <- S.thresh[corrupted.rows,corrupted.cols]
	
	def.S.thresh <- rpca.pipe.res$rob.res_def$S
	def.S.thresh[def.S.thresh < 0] <- -1	
	def.S.thresh[def.S.thresh > 0] <- 1
	def.corrupted.rows <- which(rowSums(abs(def.S.thresh)) > 0)	
	def.corrupted.cols <- which(colSums(abs(def.S.thresh[def.corrupted.rows,]))>0)
	def.corrupted.threshold <- def.S.thresh[def.corrupted.rows,def.corrupted.cols]	
	
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
	
	
	outliers <- matrix(0,nrow(X),4)
	rownames(outliers) <- rownames(X)

	#colnames(outliers) <- c("Low Rank > 90% Quant","Sparse Low Rank > 90% Quant","RPCA Corrupted")
	colnames(outliers) <- c(
	                        paste0("Low Rank > ",paste0(quantile.thresh * 100,"% Quant")),
	                        paste0("Sparse Low Rank > ",paste0(quantile.thresh * 100,"% Quant")),
	                        "Default RPCA Corrupted",
	                        "Sparse Search RPCA Corrupted")
	                        
	outliers[names(low.rank.outliers),paste0("Low Rank > ",paste0(quantile.thresh * 100,"% Quant"))] <- 1
	outliers[names(sparse.low.rank.outliers),paste0("Sparse Low Rank > ",paste0(quantile.thresh * 100,"% Quant"))] <- 1
	outliers[names(def.corrupted.rows),"Default RPCA Corrupted"] <- 1
	outliers[names(corrupted.rows),"Sparse Search RPCA Corrupted"] <- 1
	
	
	return(list(outliers=outliers,
	            corrupted.threshold= corrupted.threshold,
	            def.corrupted.threshold= def.corrupted.threshold,
		dists=list(
				X.from.S_def_OtoT_orth.dist=X.from.S_def$OtoT_orth.dist,
				L_def.from.S_def_OtoT_orth.dist=L_def.from.S_def$OtoT_orth.dist,
				X.from.S_sparse_OtoT_orth.dist=X.from.S_sparse$OtoT_orth.dist,
				L_sparse.from.S_sparse_OtoT_orth.dist=L_sparse.from.S_sparse$OtoT_orth.dist
				)
		))

	
}