tolEllipsePlot2 <- function (x, m.cov = covMcd(x), cutoff = NULL, id.n = NULL,
          tol = 1e-07, xlab = "", ylab = "", main = "Tolerance ellipse (97.5%)",
          txt.leg = c("robust", "classical"), col.leg = c("red", "blue"),
          lty.leg = c("solid", "dashed")) {

  ellips2 <- function(loc, cov) {
    dist <- sqrt(qchisq(0.975, 2))
    A <- solve(cov)
    eA <- eigen(A)
    ev <- eA$values
    lambda1 <- max(ev)
    lambda2 <- min(ev)
    eigvect <- eA$vectors[, order(ev)[2]]
    z <- seq(0, 2 * pi, by = 0.01)
    z1 <- dist/sqrt(lambda1) * cos(z)
    z2 <- dist/sqrt(lambda2) * sin(z)
    alfa <- atan(eigvect[2]/eigvect[1])
    r <- matrix(c(cos(alfa), -sin(alfa), sin(alfa), cos(alfa)),
                ncol = 2)
    elli <- t(loc + t(cbind(z1, z2) %*% r))
      ## ok this is the boundary of the ellipse.
    return(list(elli=elli,z1=z1,z2=z2,r=r,loc=loc,a1=dist/sqrt(lambda1),a2=dist/sqrt(lambda2)))
  }

  if (is.data.frame(x))
    x <- data.matrix(x)
  if (!is.matrix(x) || !is.numeric(x))
    stop("x is not a numeric dataframe or matrix.")
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (p != 2)
    stop("Dimension {= ncol(x)} must be 2!")
  if (!is.numeric(m.cov$center) || !is.numeric(m.cov$cov))
    stop("argument 'm.cov' must have numeric components 'center' and 'cov'")

  x.loc <- m.cov$center
  x.cov <- n/(n - 1) * m.cov$cov
  xM <- colMeans(x)

  e.robust.res <- ellips2(loc = xM, cov = n/(n - 1) * cov.wt(x)$cov)
  e.classic.res <- ellips2(loc = x.loc, cov = x.cov)

  z1 <- e.robust.res$elli
  z2 <- e.classic.res$elli

  x1 <- c(min(x[, 1], z1[, 1], z2[, 1]), max(x[, 1], z1[, 1],
                                             z2[, 1]))
  y1 <- c(min(x[, 2], z1[, 2], z2[, 2]), max(x[, 2], z1[, 2],
                                             z2[, 2]))

  #md <- sqrt(mahalanobis(x, xM, cov(x), tol = tol))

  rd <- sqrt(mahalanobis(x, m.cov$center, m.cov$cov, tol = tol))
  if (missing(cutoff) || is.null(cutoff))
    cutoff <- sqrt(qchisq(0.975, df = 2))
  if (missing(id.n) || is.null(id.n))
    id.n <- sum(rd > cutoff)
  plot(x, xlim = x1, ylim = y1, xlab = xlab, ylab = ylab, main = main)
  box()
  xrange <- par("usr")
  xrange <- xrange[2] - xrange[1]
  if (id.n >= 1) {
    ind <- sort(rd, index.return = TRUE)$ix[(n - id.n + 1):n]
    text(x[ind, 1] + xrange/50, x[ind, 2], ind)
  }
  points(z2, type = "l", lty = lty.leg[1], col = col.leg[1])
  #if (classic) {
    points(z1, type = "l", lty = lty.leg[2], col = col.leg[2])
    legend("bottomright", txt.leg, lty = lty.leg, col = col.leg)
  #}
  #invisible()
  return( list(ind=ind,z1=z1,z2=z2,e.robust.res=e.robust.res,e.classic.res=e.classic.res) )
}
