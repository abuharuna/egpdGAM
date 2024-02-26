
evgam_egpd <- function(formula, data, family="egpd", correctV=TRUE, rho0=0, 
                       inits=NULL, outer="fd", control=NULL, removeData=FALSE, trace=0,
                       knots=NULL, maxdata=1e20, maxspline=1e20, compact=FALSE,
                       ald.args=list(), exi.args=list(), pp.args=list(), sandwich.args=list()) {
  
## setup family
  message("****setting up data*****")
family.info <- .setup.family(family, pp.args)

if (is.null(family.info$lik.fns$d340))
  outer <- "fd"

## setup formulae
formula <- .setup.formulae(formula = formula, npar = family.info$npar, npar2 =  family.info$npar2, data, trace)
response.name <- attr(formula, "response.name")

## setup mgcv objects and data
temp.data <- .setup.data(data = data, responsename =  response.name, formula = formula, family, nms = family.info$nms, 
  removeData, exi.args, ald.args, pp.args, knots, maxdata, 
  maxspline, compact, sargs = sandwich.args, outer = tolower(outer), trace)
data <- temp.data$data

## initialise inner iteration
message("****inner iterations*****")
beta <- .setup.inner.inits(inits = inits, likdata = temp.data$lik.data, likfns = family.info$lik.fns, npar = family.info$npar, family)
lik.data <- .sandwich(temp.data$lik.data, beta)
if (trace > 0 & lik.data$adjust > 0) cat(paste("\n Sandwich correct lambda =", signif(lik.data$k, 3), "\n"))

## check whether any smoothing parameters need estimating

smooths <- length(temp.data$gotsmooth) > 0
message("****outer iterations*****")
if (smooths) {

## initialise outer iteration
S.data <- .joinSmooth(temp.data$gams)
nsp <- length(attr(S.data, "Sl"))
if (is.null(rho0)) {
    diagSl <- sapply(attr(S.data, "Sl"), diag)
    rho0 <- apply(diagSl, 2, function(y) uniroot(.guess, c(-1e2, 1e2), d=attr(beta, "diagH"), s=y)$root)
} else {
    if (length(rho0) == 1) rho0 <- rep(rho0, nsp)
}

lik.data$S <- .makeS(S.data, exp(rho0))

## perform outer iteration
fit.reml <- .outer(rho0 = rho0, beta = beta, likfns = family.info$lik.fns, likdata = lik.data, Sdata = S.data, control = control, correctV, outer = lik.data$outer, trace)

sp <- exp(fit.reml$par)
lik.data$S <- .makeS(S.data, sp)

} else {

S.data <- NULL
fit.reml <- .outer.nosmooth(beta, family.info$lik.fns, lik.data, control, trace)

}
message("****finalizing*****")
## covariance matrices
VpVc <- .VpVc(fitreml = fit.reml, likfns = family.info$lik.fns, likdata = lik.data, Sdata = S.data, correctV=correctV, sandwich=temp.data$sandwich, smooths=smooths, trace=trace)

## effective degrees of freedom
edf <- .edf(fit.reml$beta, family.info$lik.fns, lik.data, VpVc, temp.data$sandwich)

## update mgcv objects
names(temp.data$gams) <- family.info$nms
gams <- .swap(fitreml = fit.reml, gams = temp.data$gams, likdata = lik.data, VpVc, gotsmooth = temp.data$gotsmooth,edf =  edf, smooths = smooths)

## add extra things that make an evgam object
## differ from a list of mgcv objects

gams <- .finalise(gams, data, family.info$lik.fns, lik.data, S.data, fit.reml, VpVc, family, temp.data$gotsmooth, formula, response.name, removeData, edf)
message("****DONE*****")
return(gams)
}

#' @rdname evgam
#' @name fevgam
#' @export
NULL
fevgam <- function(...) {
message("`fevgam' will soon be deprecated: please migrate to `evgam'.")
evgam(...)
}

#' Extract Model Fitted Values
#'
#' @param object a fitted \code{evgam} object
#' @param ... not used
#'
#' @examples
#'
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' fitted(m_gev)
#'
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
#' 
fitted.evgam <- function(object, ...) {
predict(object)
}

#' Simulations from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a data frame
#' @param type a character string, as in \code{predict.evgam}; defaults to \code{"quantile"}
#' @param probs a scalar or vector of probabilities for quantiles; defaults to NULL
#' @param threshold a scalar, vector or matrix, which is added to each simulation if \code{family == "gpd"}; defaults to 0
#' @param marginal a logical: should simulations integrate out smoothing parameter uncertainty? Defaults to TRUE
#' @param ... arguments to be passed to \code{predict.evgam}
#'
#' @return Simulations of parameters or quantiles
#'
#' @seealso \link{predict.evgam}
#'
#' @examples
#'
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' # simulations of link GEV parameters for fremantle data
#' simulate(m_gev, nsim=5)
#' # simulations for Year 1989
#' y1989 <- data.frame(Year = 1989)
#' # link GEV parameter simulations
#' simulate(m_gev, nsim=5, newdata = y1989)
#' # GEV parameter simulations
#' simulate(m_gev, nsim=5, newdata = y1989, type = "response")
#' # 10-year return level simulations
#' simulate(m_gev, nsim=5, newdata = y1989, type= "quantile", prob = .9)
#' # 10- and 100-year return level simulations
#' simulate(m_gev, nsim=5, newdata = y1989, type= "quantile", prob = c(.9, .99))
#'
#' @export
#' 
simulate.evgam <- function(object, nsim=1e3, seed=NULL, newdata, 
  type="link", probs=NULL, threshold=0, marginal=TRUE, ...) {
if (!is.null(probs)) 
  type <- "quantile"
# if (is.null(newdata)) newdata <- object$data
if (type %in% c("link", "response")) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  family <- object$family
  if (marginal) {
    V.type <- "Vc" 
  } else {
    V.type <- "Vp"
  }
  B <- .pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
  idpars <- object$idpars
  X <- predict.evgam(object, newdata, type="lpmatrix")
  nms <- names(X)
  B <- lapply(seq_along(X), function(i) B[idpars == i, , drop=FALSE])
  X <- lapply(seq_along(X), function(i) X[[i]] %*% B[[i]])
  names(X) <- nms
  if (type == "response") {
    if (family != "exi") {
      unlink <- which(substr(nms, 1, 3) == "log")
      for (i in unlink) {
        X[[i]] <- exp(X[[i]])
        if (substr(nms[i], 1, 5) == "logit")
          X[[i]] <- X[[i]] / (1 + X[[i]])
      } 
    } else {
      X[[i]] <- object$linkfn(X[[i]])
    }
  }
nms <- gsub("cloglog", "", nms)
nms <- gsub("probit", "", nms)
nms <- gsub("logit", "", nms)
nms <- gsub("log", "", nms)
names(X) <- nms
}
if (type == "quantile") {
  X <- simulate.evgam(object, nsim, seed, newdata, "response")
  out <- list()
  for (i in seq_along(probs)) {
    if (object$family == "gpd") {
      out[[i]] <- threshold + .qgpd(probs[i], 0, X[[1]], X[[2]])
    } else {
      out[[i]] <- .qgev(probs[i], X[[1]], X[[2]], X[[3]])
    }
  }
  names(out) <- paste("q", probs, sep=":")
  if (length(probs) == 1) {
    X <- out[[1]]
  } else {
    if (nrow(X[[1]]) == 1) {
      X <- t(sapply(out, c))
    } else {
      X <- out
    }
  }
}
return(X)
}

#' Log-likelihood, AIC and BIC from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param ... not used
#'
#' @return A scalar
#'
#' @examples
#' 
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' logLik(m_gev)
#' AIC(m_gev)
#' BIC(m_gev)
#'
#' @export
#' 
logLik.evgam <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- object$logLik
attr(out, "df") <- attr(object, "df")
attr(out, "nobs") <- nobs(object)
class(out) <- "logLik"
out
}

#' Plot a fitted \code{evgam} object
#'
#' @param x a fitted \code{evgam} object
#' @param onepage logical: should all plots be on one page, or on separate pages? Defaults to \code{TRUE}
#' @param which a vector of integers identifying which smooths to plot. The default \code{NULL} plots all smooths
#' @param main a character string or vector of plot titles for each plot. If not supplied default titles are used
#' @param ask logical: ask to show next plots if too many figures for current device?
#' @param ... extra arguments to pass to \link[mgcv]{plot.gam}
#'
#' @return Plots representing all one- or two-dimensional smooths
#'
#' @examples
#'
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' plot(m_gev)
#'
#' @export
#' 
plot.evgam <- function(x, onepage = TRUE, which = NULL, main, ask = !onepage, ...) {
x <- x[x$gotsmooth]
if (is.null(which)) {
  nplot <- sum(unlist(lapply(x, function(x) as.integer(sapply(x$smooth, function(y) y$plot.me)))))
  which <- seq_len(nplot)
} else {
  nplot <- length(which)
}
if (onepage) {
  omfrow <- par("mfrow")
#   nmfrow <- n2mfrow(nplot, asp = 3) # use when R 4.0.0 is old version
  nmfrow <- rev(n2mfrow(nplot))
  par(mfrow = nmfrow)
}

if (ask) {
  oask <- par("ask")
  if (nplot > prod(par("mfrow")) && dev.interactive()) {
    par(ask = TRUE)
  } else {
    ask <- FALSE
  }
}

current <- 1
for (i in seq_along(x)) {
  for (j in seq_along(x[[i]]$smooth)) {
    if (current %in% which) {
      if (missing(main)) {
        mgcv::plot.gam(x[[i]], select = j, main = paste(names(x)[i], x[[i]]$smooth[[j]]$label, sep = ": "), ...)
      } else {
        mgcv::plot.gam(x[[i]], select = j, ...)
      }
    }
    current <- current + 1
  }
}
if (onepage)
  par(mfrow = omfrow)
if (ask)
  par(ask = oask)
}

#' Summary method for a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param ... not used
#'
#' @details
#' 
#' The key part of summary.evgam is p-values for smooths.
#' The tests use code directly taken from \code{mgcv 1.8-14}. This is 
#' to avoid use of \code{mgcv:::...} . Tests implement the method of
#' Wood (2013).
#'
#' @references
#'
#' Wood, S. N., (2013) On p-values for smooth components of an extended
#' generalized additive model, Biometrika 100(1) 221--228
#'
#' @return A \code{summary.evgam} object
#'
#' @examples
#'
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' summary(m_gev)
#'
#' @name summary.evgam
#'
#' @export
#' 
summary.evgam <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- list()
out[[1]] <- .parametric.summary.evgam(object)
out[[2]] <- .smooth.summary.evgam(object)
class(out) <- "summary.evgam"
out
}

#' @param x a \code{summary.evgam} object
#'
#' @rdname summary.evgam
#' 
#' @export
#' 
print.summary.evgam <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
cat("\n")
cat("** Parametric terms **")
tab <- lapply(x[[1]], .tidyParametricTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
cat("\n")
cat("** Smooth terms **")
tab <- lapply(x[[2]], .tidySmoothTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
invisible(x)
}

#' Print a fitted \code{evgam} object
#'
#' @param x a fitted \code{evgam} object
#' @param ... not used
#'
#' @return The call of the \code{evgam} object
#'
#' @examples
#'
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' print(m_gev)
#'
#' @export
#' 
print.evgam <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
print(x$call)
invisible(x)
}

#' Bind a list a data frames
#'
#' @param x a list of data frames
#'
#' @return A data frame
#'
#' @examples
#'
#' z <- list(data.frame(x=1, y=1), data.frame(x=2, y=2))
#' dfbind(z)
#'
#' @seealso \link[base]{rbind}
#'
#' @export
#' 
dfbind <- function(x) {
nms <- names(x[[1]])
cls <- sapply(x[[1]], class)
x <- lapply(nms, function(i) unlist(lapply(x, function(y) y[,i])))
x <- as.data.frame(x)
dt <- cls == "Date"
if (any(dt)) {
  for (i in which(dt)) x[,i] <- as.Date(x[,i], origin="1970-01-01")
}
names(x) <- nms
x
}

#' Scatter plot, with variable-based point colours
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param z a variable for defining colours
#' @param n an integer giving the number of colour levels, supplied to \link[base]{pretty}
#' @param breaks a vector or breaks for defining color intervals; defaults to \code{NULL}, so \link[base]{pretty} and \code{n} are used on \code{z}
#' @param palette a function for the color palette, or colors between \code{breaks}; defaults to \link[grDevices]{heat.colors}
#' @param rev logical: should the palette be reversed? Defaults to \code{TRUE}
#' @param pch an integer giving the plotting character, supplied to \link[graphics]{plot}
#' @param add should this be added to an existing plot? Defaults to \code{FALSE}
#' @param z.lim xxx
#' @param ... other arguments passed to \link[graphics]{plot}
#' @param legend should a legend be added? Defaults to code{FALSE}
#' @param n.legend an integer giving the approximate number of legend entries; defaults to 6
#' @param legend.pretty logical: should the legend values produced by \[base]{pretty}? Othewrwise they are exact. Defaults to \code{TRUE}
#' @param legend.x passed to \link[graphics]{legend}'s \code{x} argument
#' @param legend.y passed to \link[graphics]{legend}'s \code{y} argument
#' @param legend.horiz passed to \link[graphics]{legend}'s \code{horiz} argument
#' @param legend.bg passed to \link[graphics]{legend}'s \code{bg} argument
#' @param legend.plot passed to \link[graphics]{legend}'s \code{plot} argument
#'
#' @return A plot
#'
#' @examples
#'
#' x <- runif(50)
#' y <- runif(50)
#' colplot(x, y, x * y)
#' colplot(x, y, x * y, legend=TRUE, legend.x="bottomleft")
#' colplot(x, y, x * y, legend=TRUE, legend.pretty=FALSE, n.legend=10, 
#'   legend.x="bottomleft", legend.horiz=TRUE)
#'
#' @export
#' 
colplot <- function(x, y, z, n = 20, z.lim = NULL, breaks = NULL, palette = heat.colors, rev = TRUE, 
  pch = 21, add = FALSE, ..., legend = FALSE, n.legend = 6, legend.pretty = TRUE, 
  legend.plot = TRUE, legend.x, legend.y = NULL, legend.horiz = FALSE, legend.bg = par("bg")) {
if (!is.null(breaks)) {
  brks <- breaks
} else {
  if (!is.null(z.lim)) {
    brks <- pretty(z.lim, n)
  } else {
    brks <- pretty(z, n)
  }
}
n <- length(brks[-1])
if (inherits(palette, "function")) {
  pal <- palette(n)
  if (rev) 
    pal <- rev(pal)
} else {
  pal <- palette
  if (length(pal) != length(breaks) - 1)
    stop("length(palette) != length(breaks) - 1.")
}  
col <- pal[as.integer(cut(z, brks))]
if (!add) 
  plot(x, y, type="n", ...)
if (pch %in% 21:25) {
  points(x, y, bg=col, pch=pch, ...)
} else {
  points(x, y, col=col, pch=pch, ...)
}
if (legend) {
  if (!legend.pretty) {
    brks2 <- seq(brks[1], brks[length(brks)], l=n.legend)
  } else {
    brks2 <- pretty(brks, n.legend)
  }
  n.legend <- length(brks2[-1])
  pal2 <- palette(n.legend)
  if (rev) pal2 <- rev(pal2)
  brks2 <- format(brks2, width=3, flag="0")
  lg <- paste(brks2[1:(n.legend - 1)], brks2[2:n.legend], sep=" - ")
  if (pch %in% 21:25) {
    lg <- legend(x=legend.x, y=legend.y, legend=lg, pch=pch, pt.bg=pal2,
      horiz=legend.horiz, bg=legend.bg, plot=legend.plot)
  } else {
    lg <- legend(x=legend.x, y=legend.y, legend=lg, pch=pch, col=pal2,
      horiz=legend.horiz, bg=legend.bg, plot=legend.plot)
  }
}
if (!legend.plot) 
  lg
}

#' Moore-Penrose pseudo-inverse of a matrix
#'
#' @param x a matrix
#' @param tol a scalar
#' 
#' @details
#' 
#' This function is merely a wrapper for Armadillo's pinv function with its
#' default settings, which, in particular uses the divide-and-conquer
#' method. If \code{tol} isn't provided Armadillo's default for pinv is used.

#' \code{ginv.evgam} mimics \link[MASS]{ginv} using Armadillo's pinv.
#'
#' @return A matrix
#' 
#' @references
#'
#' http://arma.sourceforge.net/docs.html#pinv
#'
#' @seealso \link[MASS]{ginv}
#'
#' @export
#' 
pinv <- function(x, tol=-1) {
armapinv(x, tol)
}

#' @rdname pinv
#'
#' @export
#' 
ginv.evgam <- function(x, tol=sqrt(.Machine$double.eps)) {
armaginv(x, tol)
}

#' More Sequence Generation
#'
#' Generate a sequence of values between a range.
#'
#' @param x a 2-vector
#' @param length an integer
#'
#' @return A vector
#'
#' @seealso \link[base]{seq}, \link[base]{seq_len}, \link[base]{seq_along}
#'
#' @examples
#'
#' seq_between(c(1, 9))
#' seq_between(range(runif(10)), 5)
#'
#' @export
#' 
seq_between <- function(x, length=NULL) {
if (is.null(length)) {
    return(seq(x[1], x[2]))
    } else {
    return(seq(x[1], x[2], length=length))
    }
}
