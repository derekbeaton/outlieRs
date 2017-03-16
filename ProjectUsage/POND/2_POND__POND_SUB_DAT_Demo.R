require(outlieRs)

## These are comments
## Any lines preceded by a # will be ignored by R, so you can write whatever you'd like!

############################################################################
## Neuroimaging example ##
############################################################################

## We recommend for now that any preprocessing, scaling, and correction factors be done *before hand*. Read in a "cleaned" data file for outlier analyses
## The cleaning/correcting/scaling processes should be done in whatever software you are comfortable with.

# Read data in from .csv file
load('POND.SUB.DAT.rda')
load('POND.SUB.DEMOS.rda')

data.for.outliers <- expo.scale(POND.SUB.DAT)



# What are the dimensions of the data?
dim(data.for.outliers)
  ## specifically get rows and columns
  nrow(data.for.outliers)
  ncol(data.for.outliers)

# What are the column names?
  colnames(data.for.outliers)
# What are the row names?
  rownames(data.for.outliers)

# What do the first three rows and columns look like?
data.for.outliers[1:3,1:3]

# What do the first three rows and *all* the columns look like?
data.for.outliers[1:3,]

# Let's get some descriptors and summaries of the data.
str(data.for.outliers)
summary(data.for.outliers)



##### PCA Diagnostics #####
pca.diagnosticsResult <- pca.diagnostics(data.for.outliers,center = F,scale=F)



  ### aggregate outlier results below.

##### Minimum covariance determinant (MCD) #####
# to see the help file (e.g., short description, inputs, outputs), type the following line (without the #) into the console:
# ?multiOut
mcdResult <- multiOut(data.for.outliers,rmdo_alpha=0.9)

# looking at results
## this provides a full set of results
names(mcdResult)
## this is the binary decision: if a row is outlier or not
mcdOut <- mcdResult$outlier_decision
## look at the first few
head(mcdOut)

## diagnostic plots
#printOuts(mcdResult)


##### MCD + Bootstrap PCA (for Mahalanobis and score distances; Boot-MCD) #####
# looking at results
## this provides a full set of results
bootResult <- mcd.boot_find.outliers(data.for.outliers)
## this is the binary decision: if a row is outlier or not
bootOut <- bootResult$outlier_decision
## look at the first few
head(bootOut)


##### Candes' Robust PCA #####
# # does not work for this particular example; ill suited for n >> p
# # looking at results
#   ## this provides a full set of results
# rpcaResult <- rpca_find.outliers(data.for.outliers)
#   ## this is the binary decision: if a row is outlier or not
# rpcaOut <- rpcaResult$outlier_decision
#   ## look at the first few
# head(rpcaOut)


# combine results across all three methods
## Get the rownames of individuals.
these.rows <- rownames(mcdOut)
## Use cbind to column bind vectors. Index them by 'these.rows'
allOut <- cbind(
  mcdOut[these.rows,'mcd_outlier'],
  bootOut[these.rows,c('mah_outlier','chi_outlier')]
)
head(allOut)
colnames(allOut) <- c("MCD","BootMCD - MD","BootMCD - ChiD")

# counting the number of thresholds that identified each observation an outlier
allOutSums <- rowSums(allOut)
head(allOutSums)
# frequency of the number of thresholds
table(allOutSums)
# frequency of observations that are outliers for each threshold
colSums(allOut)
# a comparison of the number of outliers both MCD and RPCA detected vs one or the other
crossprod(as.matrix(allOut))


