require(outlieRs)

## These are comments
## Any lines preceded by a # will be ignored by R, so you can write whatever you'd like!

############################################################################
## Deidentified example ##
############################################################################

## We recommend for now that any preprocessing, scaling, and correction factors be done *before hand*. Read in a "cleaned" data file for outlier analyses
## The cleaning/correcting/scaling processes should be done in whatever software you are comfortable with.

load('fin.scaled.dat.rda')


# What are the dimensions of the data?
dim(fin.scaled.dat)
  ## specifically get rows and columns
  nrow(fin.scaled.dat)
  ncol(fin.scaled.dat)

# What are the column names?
  colnames(fin.scaled.dat)
# What are the row names?
  rownames(fin.scaled.dat)

# What do the first three rows and columns look like?
fin.scaled.dat[1:3,1:3]

# What do the first three rows and *all* the columns look like?
fin.scaled.dat[1:3,]

# Let's get some descriptors and summaries of the data.
str(fin.scaled.dat)
summary(fin.scaled.dat)




##### Minimum covariance determinant (MCD) #####
# to see the help file (e.g., short description, inputs, outputs), type the following line (without the #) into the console:
# ?multiOut
mcdResult <- multiOut(fin.scaled.dat,rmdo_alpha=0.8)

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
bootResult <- mcd.boot_find.outliers(fin.scaled.dat)
  ## this is the binary decision: if a row is outlier or not
bootOut <- bootResult$outlier_decision
  ## look at the first few
head(bootOut)


##### Candes' Robust PCA #####
# looking at results
  ## this provides a full set of results
rpcaResult <- rpca_find.outliers(fin.scaled.dat)
  ## this is the binary decision: if a row is outlier or not
rpcaOut <- rpcaResult$outlier_decision
  ## look at the first few
head(rpcaOut)


# combine results across all three methods
  ## Get the rownames of individuals.
#these.rows <- rownames(mcdOut)
these.rows <- rownames(rpcaOut)
  ## Use cbind to column bind vectors. Index them by 'these.rows'
allOut <- cbind(
                #mcdOut[these.rows,'mcd_outlier'],
                bootOut[these.rows,c('mah_outlier','chi_outlier')],
                rpcaOut[these.rows,'rpca_outlier']
              )
head(allOut)
colnames(allOut) <- c("MCD","BootMCD - MD","BootMCD - ChiD","RPCA")

# counting the number of thresholds that identified each observation an outlier
allOutSums <- rowSums(allOut)
head(allOutSums)
# frequency of the number of thresholds
table(allOutSums)
# frequency of observations that are outliers for each threshold
colSums(allOut)
# a comparison of the number of outliers both MCD and RPCA detected vs one or the other
crossprod(as.matrix(allOut))

