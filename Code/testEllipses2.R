ellips2 <- function(loc, cov) {
  dista <- sqrt(qchisq(0.975, 2))
  A <- solve(cov)
  eA <- eigen(A)
  ev <- eA$values
  lambda1 <- max(ev)
  lambda2 <- min(ev)
  eigvect <- eA$vectors[, order(ev)[2]]
  z <- seq(0, 2 * pi, by = 0.01)
  z1 <- dista/sqrt(lambda1) * cos(z)
  z2 <- dista/sqrt(lambda2) * sin(z)
  alfa <- atan(eigvect[2]/eigvect[1])
  r <- matrix(c(cos(alfa), -sin(alfa), sin(alfa), cos(alfa)),ncol = 2)
  elli <- t(loc + t(cbind(z1, z2) %*% r))
  ## ok this is the boundary of the ellipse.
  return(list(elli=elli,z1=z1,z2=z2,r=r,loc=loc,a1=dista/sqrt(lambda1),a2=dista/sqrt(lambda2)))
}



data(hbk)
dat <- data.matrix(hbk[, 1:2])
mcd <- covMcd(dat)       # compute mcd in advance

mcd.center <- mcd$center
mcd.cov <- mcd$cov

data.center <- colMeans(dat)
data.cov <- cov(dat)


#e.robust.res <- ellips2(loc = mcd.center, cov = nrow(dat)/(nrow(dat) - 1) * cov.wt(dat)$cov)
e.robust2.res <- ellips2(loc = mcd.center, cov = mcd.cov)
e.classic.res <- ellips2(loc = data.center, cov = data.cov)


Z1 <- pointsToEllipsoid(dat,mcd.cov,mcd.center)
Z1.test <- ellipseInOut(Z1,p=.95)

#z1 <- e.robust.res$elli
z1 <- e.robust2.res$elli
z2 <- e.classic.res$elli

x1 <- c(min(dat[, 1], z1[, 1], z2[, 1]), max(dat[, 1], z1[, 1], z2[, 1]))
y1 <- c(min(dat[, 2], z1[, 2], z2[, 2]), max(dat[, 2], z1[, 2], z2[, 2]))


plot(dat, xlim = x1, ylim = y1)
text(dat,labels=1:nrow(dat),pos=1)
#points(z1, type = "l", lty = 1, col="red")
points(z1, type = "l", lty = 1, col="green")
points(z2, type = "l", lty = 1, col="blue")




((dat[14,1] - data.center[1])^2 / (e.classic.res$a1^2) + (dat[14,2] - data.center[2])^2 / (e.classic.res$a2^2))


for(i in 1:nrow(dat)){

  print(  ((dat[i,1] - colMeans(e.classic.res$elli)[1])^2 / (e.classic.res$a1^2)) + ((dat[i,2] - colMeans(e.classic.res$elli)[2])^2 / (e.classic.res$a2^2)) )

}


for(i in 1:nrow(dat)){

  print( unique(rowSums(abs(dat[i,]) > abs(e.classic.res$elli))) )

}



