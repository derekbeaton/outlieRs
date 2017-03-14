rm(list=ls())

data(hbk)
dat <- scale(data.matrix(hbk[, 1:2]),scale=F)


mcd <- covMcd(dat)       # compute mcd in advance

mcd.center <- mcd$center
mcd.cov <- mcd$cov

data.center <- colMeans(dat)
data.cov <- cov(dat)

#e.robust.res <- ellips2(loc = mcd.center, cov = mcd.cov)
#e.classic.res <- ellips2(loc = data.center, cov = data.cov)


#plot(dat[,2]~dat[,1])


z1 <- pointsToEllipsoid(dat,mcd.cov,mcd.center)
z1.test <- ellipseInOut(z1,p=.95)

z2 <- pointsToEllipsoid(dat,data.cov,data.center)
z2.test <- ellipseInOut(z2,p=.95)


x1 <- c(-max(abs(dat[,1])),max(abs(dat[,1])))
y1 <- c(-max(abs(dat[,2])),max(abs(dat[,2])))

  ### this is the solution.
plot(dat, xlim = x1, ylim = y1,pch=20,col="grey80")
tmp <- addEllipse(mcd.center,mcd.cov,p.interval = .95,col="red",lty=1)
tmp2 <- addEllipse(data.center,data.cov,p.interval = .95,col="blue",lty=2)
points(dat[!z1.test,],bg="red",pch=21)
  text(dat[!z1.test,],labels=rownames(dat[!z2.test,]),pos=1,col="red")
points(dat[!z2.test,],bg="blue",pch=21)
  text(dat[!z2.test,],labels=rownames(dat[!z2.test,]),pos=1,col="blue")


library(ExPosition)
source('R/pca.diagnostics.R')
res <- pca.diagnostics(scale(data.matrix(hbk[, 1:2]),scale=F),center = F,scale = F)

pca.res <- epPCA(scale(data.matrix(hbk[, 1:2]),scale=F),F,F,graphs=F)

