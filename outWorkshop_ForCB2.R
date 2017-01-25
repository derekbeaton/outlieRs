require(outlieRs)

############################################################################
## Neuroimaging example ##
############################################################################

# Importing the data
nidat <- read.csv("Clean_FS_a_data.csv",header=T)
dim(nidat)
nidat[1:3,1:3]

colnames(nidat)

str(nidat)
summary(nidat)


# Data Prep
# Matrix must be numbers only

# moving identifier into rownames
nidat1 <- nidat[,-1] # removing first column
rownames(nidat1) <- nidat$Measure.volume

nidat1 <- as.matrix(nidat1)
summary(nidat1)
nidat1[1:3,1:3]


# For ONDRI:

# centered and scaled each variable
# imputing missing values with the column mean ---> data must be complete
# considered proportion of missingness per observation and variable
### if a variable was missing 90% of values, we excluded it
nidat2 <- apply(nidat1,2,function(i){(i-mean(i,na.rm=T))/sd(i,na.rm=T)}) 


##### MCD! #####
?multiOut
mcdResult <- multiOut(nidat2[,1:14],rmdo_alpha=0.9)

# looking at results
names(mcdResult)
mcdOut <- mcdResult$outlier_decision
head(mcdOut)

printOuts(mcdResult)

##### Boot-MCD! #####

bootResult <- mcd.boot_find.outliers(nidat2)

bootOut <- bootResult$outlier_decision
head(bootOut)

# RPCA
rpcaResult <- rpca_find.outliers(nidat2)

rpcaOut <- rpcaResult$outlier_decision
head(rpcaOut)


# looking at results across all three methods
# merge function combines tables based on common column names (why I produced data frames with an ID column)
allOut0 <- merge(mcdOut,bootOut)
allOut <- merge(allOut0,rpcaOut)
head(allOut)

allOut1 <- allOut[,-1]
rownames(allOut1) <- allOut$ID

# counting the number of thresholds that identified each observation an outlier
allOutSums <- rowSums(allOut1)
head(allOutSums)
# frequency of the number of thresholds
table(allOutSums)
# propotion of observations that are outliers for each threshold
colMeans(allOut1)
# a comparison of the number of outliers both MCD and RPCA detected vs one or the other
table(allOut1$mcd_outlier,allOut1$rpca_outlier)

