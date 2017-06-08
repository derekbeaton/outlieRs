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


##### PCA Diagnostics #####
pca.diagnosticsResult <- pca.diagnostics(fin.scaled.dat,center = F,scale=F)

  ## robust ellipse outliers
which(pca.diagnosticsResult$ellipse.outliers$robust.outs)
  ## classic ellipse outliers
which(pca.diagnosticsResult$ellipse.outliers$classic.outs)

  ### quantile-based outliers
  ## Chi2 Distance
pca.diagnosticsResult$quantile.outliers$out.quantile.cd
  ## Mahalanobis Distance
pca.diagnosticsResult$quantile.outliers$out.quantile.md
  ## Chi2 & Mahalanobis Distance
pca.diagnosticsResult$quantile.outliers$out.quantile.both






