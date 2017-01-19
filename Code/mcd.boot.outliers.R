mcd.boot.outliers <- function(dist.mat,full.cloud.thresh=.99,frequency.thresh=NA,dims=NA,mcd.thresh=NA){

		
		num.samps <- ncol(dist.mat) 
		
		if(is.na(frequency.thresh)){
			frequency.thresh <- 1-full.cloud.thresh
		}
		do.mcd <- T
		if(!is.na(dims)){
			if(is.na(mcd.thresh)){
				mcd.thresh <- .75	##I don't know why.
			}
		}else{
			do.mcd <- F
		}
		
		dist.vec <- c(dist.mat)
		frequency.outside <- rowSums(dist.mat > sort(dist.vec)[round(c(full.cloud.thresh) * length(dist.vec))])
		frequency.outsiders <- ifelse(frequency.outside > (num.samps * frequency.thresh),T,F)
	
		if(do.mcd){
			mass.quantile <- min(h.alpha.n(mcd.thresh,dims[1],dims[2]),dims[1]-1)
			mcd.sample.order <- sort(order(frequency.outside)[1:mass.quantile])			## in case there are no labels.
			return(list(freq.boot.outliers= frequency.outside/num.samps, thresh.boot.outliers=frequency.outsiders,mcd.style.best=mcd.sample.order))
		}
	
		return(list(
					freq.boot.outliers = frequency.outside/num.samps, 
					thresh.boot.outliers = frequency.outsiders
				))
		
	}	
	