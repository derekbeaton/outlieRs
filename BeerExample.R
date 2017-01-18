##script should be:

	##load in libraries, functions, and data required.
	
## Read in NP and Clin data sets so I can build a master set of data.
rm(list=ls())
gc()

library(rpca)
library(TExPosition)
library(MASS)
source('fast.pc.dists.R')
source('rpca.pipeline.R')
source('rpca.outliers.R')
source('mcd.boot.find.outliers.R')
source('mcd.boot.pipeline.R')
source('mcd.boot.outliers.R')

data(beer.tasting.notes)

	##we should center/scale here now.
the.data <- expo.scale(beer.tasting.notes$data)

the.data_corrected_for_ABV <- apply(the.data,2,
										function(x){
											resid(lm(x~beer.tasting.notes$sup.data[,"ABV"]))
										}
									)


### both of these functions just return the core information needed to estimate outliers.

## MCD -- these return a matrix.
beer_MCD.boot_full.combined_95 <- mcd.boot.find.outliers(the.data,center=F,scale=F)
beer.corrected_MCD.boot_full.combined_95 <- mcd.boot.find.outliers(the.data_corrected_for_ABV,center=F,scale=F)

	##for MCD: there are Chi2 and Mahalanobis distances, percentiles, and then a threshold of outlier/not outlier.
		

## RPCA -- these return lists (see below for the important stuff)
beer_rpca <- rpca.outliers(the.data)
beer.corrected_rpca <- rpca.outliers(the.data_corrected_for_ABV)

	##for RPCA: $outliers tells us who the outliers are based on quantiles
	##for RPCA: $corrupted.threshold which rows are corrupted at which columns and which direction