require(outlieRs)

## These are comments
## Any lines preceded by a # will be ignored by R, so you can write whatever you'd like!

############################################################################
## Neuroimaging example ##
############################################################################

## We recommend for now that any preprocessing, scaling, and correction factors be done *before hand*. Read in a "cleaned" data file for outlier analyses
## The cleaning/correcting/scaling processes should be done in whatever software you are comfortable with.

# Read data in from .csv file
ni.demographics <- read.csv("demographics.csv",header=T,row.names = 1)
ni.data <- read.csv("simulated.brain_z.csv",header=T,row.names = 1)


# What are the dimensions of the data?
dim(ni.data)
  ## specifically get rows and columns
  nrow(ni.data)
  ncol(ni.data)

# What are the column names?
  colnames(ni.data)
# What are the row names?
  rownames(ni.data)

# What do the first three rows and columns look like?
ni.data[1:3,1:3]

# What do the first three rows and *all* the columns look like?
ni.data[1:3,]

# Let's get some descriptors and summaries of the data.
str(ni.data)
summary(ni.data)


##### PCA Diagnostics #####
pca.diagnosticsResult <- pca.diagnostics(ni.data,center = F,scale=F)

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






