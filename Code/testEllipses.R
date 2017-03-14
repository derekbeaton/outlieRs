
source('R/tolEllipsePlot2.R')

data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
mcd <- covMcd(hbk.x)       # compute mcd in advance
## must be a 2-dimensional data set: take the first two columns :
test.rs <- tolEllipsePlot2(hbk.x[,1:2])




## an "impressive" example:
# data(telef)
# tolEllipsePlot(telef)
