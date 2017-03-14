#' Quick PCA diagnostics for outlier detection.
#'
#' @export
#'
#' @param X Data matrix
#' @param center a boolean, vector, or string. See ?expo.scale for details.
#' @param scale a boolean, vector, or string. See ?expo.scale for details.
#' @param ellipse.alpha the alpha level for classic and robust ellipses (around MD and CD)
#' @param quantile.alpha the alpha level for quantile cutoffs (for MD an CD)
#' @param graphs a boolean. If TRUE, two graphs will be shown.
#'
#' @return
#'  some lists.
#'
#' @examples
#'  ## a very stupid example for now.
#'  data(hbk)
#'  res <- pca.diagnostics(scale(data.matrix(hbk[, 1:2]),scale=F),center = F,scale = F)
#'

pca.diagnostics <- function(X, center=T, scale=T, ellipse.alpha=.95, quantile.alpha=.75, graphs=T){

  pca.res <- epPCA(X,center=center,scale=scale,graphs=F)

  chi.d <- rowSums(pca.res$ExPosition.Data$fi^2)
    max.abs.cd <- max(abs(chi.d))
  mahal.d <- rowSums(pca.res$ExPosition.Data$pdq$p^2)
    max.abs.md <- max(abs(mahal.d))

  dat <- cbind(chi.d,mahal.d)
    ## to get an estimate of just CD v MD robustness.
  mcd <- covMcd(dat)
  mcd.center <- mcd$center
  mcd.cov <- mcd$cov
  data.center <- colMeans(dat)
  data.cov <- cov(dat)

  quantile.cols <- rep("grey80",nrow(dat))
  quantile.sizes <- rep(.5,nrow(dat))
  quantile.dots <- rep(20,nrow(dat))
  cd.quant <- quantile(chi.d,probs=quantile.alpha)
  out.quantile.cd <- which(chi.d > cd.quant)
    if(length(out.quantile.cd) > 0){
      quantile.cols[out.quantile.cd] <- "mediumorchid4"
      quantile.sizes[out.quantile.cd] <- 1
      quantile.dots[out.quantile.cd] <- 21
    }
  md.quant <- quantile(mahal.d,probs=quantile.alpha)
  out.quantile.md <- which(mahal.d > md.quant)
    if(length(out.quantile.md) > 0){
      quantile.cols[out.quantile.md] <- "olivedrab3"
      quantile.sizes[out.quantile.md] <- 1
      quantile.dots[out.quantile.md] <- 21
    }
  out.quantile.both <- intersect(out.quantile.cd,out.quantile.md)
    if(length(out.quantile.both) > 0){
      quantile.cols[out.quantile.both] <- "firebrick3"
      quantile.sizes[out.quantile.both] <- 1
      quantile.dots[out.quantile.both] <- 21
    }

  #print(out.quantile.cd)
  #print(out.quantile.md)
  #print(out.quantile.both)


  ## also the list of outliers -- this will help out massivel for ellipses.
  z1 <- pointsToEllipsoid(dat,mcd.cov,mcd.center)
  z1.outs <- ellipseInOut(z1,p=ellipse.alpha)

  z2 <- pointsToEllipsoid(dat,data.cov,data.center)
  z2.outs <- ellipseInOut(z2,p=ellipse.alpha)

  ## three (or four) pictures
    ## lims
  x1 <- c(-max.abs.cd*.05,max.abs.cd)
  y1 <- c(-max.abs.md*.05,max.abs.md)

  if(graphs){
    dev.new()
    plot(dat, xlim = x1, ylim = y1,pch=20,col="grey80", main=paste0("Outside of ellipse.alpha = ",ellipse.alpha), xlab="Chi-squared Distances", ylab="Mahalanobis Distances",cex=.5)
    tmp <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="red",lty=1)
    tmp2 <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="blue",lty=2)
    points(dat[!z1.outs,],bg="red",pch=21,cex=1)
      text(dat[!z1.outs,],labels=rownames(dat[!z2.outs,]),pos=1,col="red")
    points(dat[!z2.outs,],bg="blue",pch=21,cex=1)
      text(dat[!z2.outs,],labels=rownames(dat[!z2.outs,]),pos=1,col="blue")
    legend("bottomright",legend=c("Classic ellipse","Robust ellipse"), col=c("blue","red"), lty=c(2,1))

    dev.new()
    plot(dat, xlim = x1, ylim = y1, main=paste0("Outside of quantile.alpha = ",quantile.alpha), pch=20, col="white", xlab="Chi-squared Distances", ylab="Mahalanobis Distances")
    abline(v=cd.quant,col="mediumorchid3",lwd=2,lty=2)
    abline(h=md.quant,col="olivedrab3",lwd=2,lty=2)
    points(dat,bg=quantile.cols, cex=quantile.sizes, pch=21)
  }

  return(list(ellipse.outliers=list(robust.outs=z1.outs,classic.outs=z2.outs),quantile.outliers=list(out.quantile.cd=out.quantile.cd,out.quantile.md=out.quantile.cd,out.quantile.both=out.quantile.both)))

}
