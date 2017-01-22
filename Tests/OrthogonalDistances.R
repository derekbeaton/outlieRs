rm(list=ls())
gc()
library(ExPosition)
library(chemometrics)
library(rrcov)


## use the data I am most familiar with.
data(beer.tasting.notes)
dat <- expo.scale(beer.tasting.notes$data)

num.comps <- 5

## ExPo PCA
e.res <- epPCA(dat,F,F,graphs=F)
  e.loadings <- e.res$ExPosition.Data$pdq$q[,1: num.comps]
  dat.rebuild <- e.res$ExPosition.Data$pdq$p[,1: num.comps] %*% diag(e.res$ExPosition.Data$pdq$Dv[1: num.comps]) %*% t(e.loadings)
  
## base PCA
p.res <- princomp(dat)
  p.loadings <- p.res$loadings[,1: num.comps]


## SVD PCA
svd.res <- svd(dat)
loadings<-svd.res$v[,1: num.comps]
recon1 <- dat%*%loadings%*%t(loadings)


### Now get distances.  

## rrcov
rrcov.pca <- PcaClassic(dat,k=num.comps)
rrcov.dists <- pca.distances(rrcov.pca, dat, rankMM(dat))


## 'chemometrics'
chemo.res <- pcaDiagplot(dat, p.res, a = num.comps, scale = F, plot=F,center=F)

## mine
e.dist <- rowSums((dat - dat.rebuild)^2)

## stack exchange: http://stats.stackexchange.com/questions/69304/how-is-orthogonal-distance-computed
s.dist <-rowSums( (dat-recon1)^2 )
#orthDist<-sqrt(rowSums(orthDist*orthDist))


	##OK I know why. the PcaHubert is a ROBpca object. Silly.
		### got it -- change to PCA classic.
rrcov.dists@od / chemo.res$ODist
rrcov.dists@od / sqrt(e.dist)
rrcov.dists@od / sqrt(s.dist)


chemo.res$ODist / sqrt(e.dist)
chemo.res$ODist / sqrt(s.dist)

e.dist / s.dist


