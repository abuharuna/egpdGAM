## Generalised Pareto negative log-likelihood functions

.egpd.d0 <- function(pars, likdata) {
  if (!likdata$censored) {
    out <- egpdd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
  } else {
    id <- likdata$cens.id
     out1 <- egpdd0(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 2], likdata$dupid[!id], likdata$duplicate)
     out2 <- egpdcd0(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id,1 , drop=FALSE], likdata$dupid[id], likdata$duplicate)
     out <- out1 + out2
  }
  out
}

.egpd.d12 <- function(pars, likdata) {
  if (!likdata$censored) {
    out <- egpdd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
  } else {
    ## old for interval 
    id <- likdata$cens.id
    out <- matrix(0, likdata$nobs, 9)
    out[!id, ] <- egpdd12(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 2], likdata$dupid[!id], likdata$duplicate)
    out[id, ] <- egpdcd12(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id,1 , drop=FALSE], likdata$dupid[id], likdata$duplicate)
    #******** NEw  for left censoring
    # id <- likdata$cens.id
    # outa <- outb <- matrix(0, likdata$nobs, 9)
    # outa[, ] <- egpdd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,2], likdata$dupid, likdata$duplicate)
    # outb[, ] <- egpdcd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
    # 
    # out = outa +outb
  }
  out
}

.egpd.d34 <- function(pars, likdata) {
  if (!likdata$censored) {
    out <- egpdd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
  } else {
    #**** old for interval 
    id <- likdata$cens.id
    out <- matrix(0, likdata$nobs, 25)
    out[!id, ] <- egpdd34(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 2], likdata$dupid[!id], likdata$duplicate)
    out[id, ] <- egpdcd34(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id,1 , drop=FALSE], likdata$dupid[id], likdata$duplicate)

    ##******** new for left censoring 
    # id <- likdata$cens.id
    # outa <- outb <- matrix(0, likdata$nobs, 25)
    # outa[, ] <- egpdd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,2], likdata$dupid, likdata$duplicate)
    # outb[, ] <- egpdcd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
    # 
    # out = outa +outb
  }
  out
}

.egpdfns <- list(d0=.egpd.d0, d120=.egpd.d12, d340=NULL)
####################### modify this for eGPD
# .qgpd <- function(p, loc, scale, shape, zeta=1, theta=1, m=1) {
# shape <- sign(shape) * pmax(abs(shape), 1e-6)
# out <- 1 - p ^ (1 / (m * theta))
# loc + scale * ((out / zeta)^(-shape) - 1) / shape
# }
# 
# .pgpd <- function(x, loc, scale, shape, tau, NAOK=FALSE, log=FALSE) {
# # function to evaluate daily cdf, e.g. Coles (2001, pp.138)
# below <- x < loc
# x <- pmax(x, loc)
# temp <- 1 + shape * (x - loc) / scale
# if (!NAOK) temp <- pmax(temp, 0)
# out <- temp ^ (-1/shape)
# out <- 1 - (1 - tau) * out
# if (log) out <- log(out)
# if (NAOK) out[below] <- NA
# out
# }
# 
# .dqgpd <- function(p, lscale, shape) {
# shape <- sign(shape) * pmax(abs(shape), 1e-6)
# .e1 <- 1 - p
# .e2 <- .e1^shape
# .e4 <- 1/.e2 - 1
# .e5 <- exp(lscale)
# d1 <- .e4 * .e5/shape
# d2 <- -((.e4/shape + log(.e1)/.e2) * .e5/shape)
# cbind(d1, d2)
# }
