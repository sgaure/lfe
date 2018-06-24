# $Id: bccorr.R 2004 2016-04-21 10:31:20Z sgaure $
narowsum <- function(x, group) {
  opt <- options(warn=-1)
  res <- try(rowsum(x,group), silent=TRUE)
  options(opt)
  narow <- is.na(rownames(res))
  if(any(narow))
      res[!narow,,drop=FALSE]
  else
      res
}

# try without reference
narowsum <- rowsum
# to estimate the bias of the D-variance we need to find tr(D' P_1 D (D' P_{X,F} D)^{-1})
# where P_1 projects out the mean, P_{X,F} projects out X and F. Let's make function
# which evaluates the matrix inside the trace on a vector.

# Now, note that for matrix x, D %*% x is the same as x[f1,]
# whereas t(D) %*% x is the same as rowsum(x,f1). Similarly with F and f2.

# bias corrected variances

# Note that this file uses the scalecols() function which scales
# its argument in-place, and returns it.  This is a violation of Rs
# immutable semantics. Use with caution.  a <- scalecols(b,r) will scale a, but also change b!


#' Compute limited mobility bias corrected correlation between fixed effects
#' @concept Limited Mobility Bias
#'
#' @description
#' With a model like \eqn{y = X\beta + D\theta + F\psi + \epsilon}, where \eqn{D}
#' and \eqn{F} are matrices with dummy encoded factors, one application of \pkg{lfe} is to
#' study the correlation \eqn{cor(D\theta, F\psi)}.  However, if we use
#' estimates for \eqn{\theta} and \eqn{\psi}, the resulting correlation is biased.
#' The function \code{bccorr} computes a bias corrected correlation
#' as described in \cite{Gaure (2014)}.
#' @param est an object of class '"felm"', the result of a call to
#'  \code{\link{felm}(keepX=TRUE)}.
#' @param alpha a data frame, the result of a call to \code{\link{getfe}}.
#' @param corrfactors integer or character vector of length 2. The factors to
#'    correlate. The default is fine if there are only two factors in the model.
#' @param nocovar logical. Assume no other covariates than the two
#'  factors are present, or that they are uncorrelated with them.
#' @param tol The absolute tolerance for the bias-corrected correlation.
#' @param maxsamples Maximum number of samples for the trace sample means estimates
#' @param lhs character. Name of left hand side if multiple left hand sides.

#' @return
#'  \code{bccorr} returns a named integer vector with the following fields:
#'
#'  \item{corr}{the bias corrected correlation.}
#'  \item{v1}{the bias corrected variance for the first factor specified
#'  by \code{corrfactors}.}
#'  \item{v2}{the bias corrected variance for the second factor.}
#'  \item{cov}{the bias corrected covariance between the two factors.}
#'  \item{d1}{the bias correction for the first factor.}
#'  \item{d2}{the bias correction for the second factor.}
#'  \item{d12}{the bias correction for covariance.}
#'
#'  The bias corrections have been subtracted from the bias estimates.
#'  E.g. v2 = v2' - d2, where v2' is the biased variance.

#' @details
#' The bias expressions from \cite{Andrews et al.} are of the form \eqn{tr(AB^{-1}C)}
#' where \eqn{A}, \eqn{B}, and \eqn{C} are  matrices too large to be handled
#' directly. \code{bccorr} estimates the trace by using the formula \eqn{tr(M) = E(x^t M x)}
#' where x is a vector with coordinates drawn uniformly from the set \eqn{\{-1,1\}}.
#' More specifically, the expectation is estimated by
#' sample means, i.e. in each sample a vector x is drawn, the
#' equation \eqn{Bv = Cx} is solved by a conjugate gradient method, and the
#' real number \eqn{x^t Av} is computed. 
#' 
#' There are three bias corrections, for the variances of \eqn{D\theta} (\code{vD}) and
#' \eqn{F\psi} (\code{vF}), and their covariance (\code{vDF}).The correlation is computed as
#' \code{rho <- vDF/sqrt(vD*vF)}.  The variances are estimated to a
#' relative tolerance specified by the argument \code{tol}. The covariance
#' bias is estimated to an absolute tolerance in the correlation \code{rho}
#' (conditional on the already bias corrected \code{vD} and \code{vF}) specified by
#' \code{tol}.  The CG algortithm does not need to be exceedingly precise,
#' it is terminated when the solution reaches a precision which is
#' sufficient for the chosen precision in \code{vD, vF, vDF}.
#' 
#' If \code{est} is the result of a weighted \code{\link{felm}} estimation,
#' the variances and correlations are weighted too.

#' @note
#' Bias correction for IV-estimates are not supported as of now.
#' 
#' Note that if \code{est} is the result of a call to \code{\link{felm}}
#' with \code{keepX=FALSE} (the default), the correlation will be computed
#' as if the covariates X are independent of the two factors. This will be
#' faster (typically by a factor of approx. 4), and possibly wronger.
#' 
#' Note also that the computations performed by this function are
#' non-trivial, they may take quite some time.  It would be wise to start
#' out with quite liberal tolerances, e.g. \cite{tol=0.1}, to
#' get an idea of the time requirements.
#' 
#' The algorithm used is not very well suited for small datasets with only
#' a few thousand levels in the factors.

#' @seealso \code{\link{fevcov}}

#' @examples
#' x <- rnorm(500)
#' x2 <- rnorm(length(x))
#' 
#' ## create individual and firm
#' id <- factor(sample(40,length(x),replace=TRUE))
#' firm <- factor(sample(30,length(x),replace=TRUE,prob=c(2,rep(1,29))))
#' foo <- factor(sample(20,length(x),replace=TRUE))
#' ## effects
#' id.eff <- rnorm(nlevels(id))
#' firm.eff <- rnorm(nlevels(firm))
#' foo.eff <- rnorm(nlevels(foo))
#' ## left hand side
#' y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + foo.eff[foo] + rnorm(length(x))
#' 
#' # make a data frame
#' fr <- data.frame(y,x,x2,id,firm,foo)
#' ## estimate and print result
#' est <- felm(y ~ x+x2|id+firm+foo, data=fr, keepX=TRUE)
#' # find bias corrections
#' bccorr(est)

#' @references
#'  Gaure, S. (2014), \cite{Correlation bias correction in two-way
#'  fixed-effects linear regression}, Stat 3(1):379:390, 2014.
#' @export
bccorr <- function(est, alpha=getfe(est), corrfactors=1L:2L,
                   nocovar=(length(est$X)==0) && length(est$fe)==2,
                   tol=0.01, maxsamples=Inf, lhs=NULL) {

  if(nlevels(est$cfactor) > 1) stop('Data should have just a single connected component')
  if(length(est$fe) == 2 && nlevels(est$cfactor) != 1) stop('Bias correction only makes sense on data with 1 component')
  if(length(est$fe) < 2) stop('Bias correction only makes sense for two factors')
  if(length(corrfactors) != 2) stop('corrfactors must have length 2')
  if(is.character(corrfactors)) corrfactors <- match(corrfactors, names(est$fe))
  if(min(corrfactors) < 1 || max(corrfactors) > length(est$fe))
      stop('corrfactors specifies too small or large index')
  if(!('fe' %in% colnames(alpha))) stop('alpha must have an "fe" column')
#  if(!is.null(est$weights)) warning("Bias correction with weights not yet fully correct")  

#  fe <- est$fe
  f1 <- est$fe[[corrfactors[[1]]]]
  f2 <- est$fe[[corrfactors[[2]]]]
  nf1 <- names(est$fe)[[corrfactors[[1]]]]
  nf2 <- names(est$fe)[[corrfactors[[2]]]]
  if(!is.null(attr(f1,'x',exact=TRUE))) stop('Interacted factors "',nf1,'" are not supported')
  if(!is.null(attr(f2,'x',exact=TRUE))) stop('Interacted factors "',nf2,'" are not supported')
  effnam <- 'effect'
  if(length(est$lhs) == 1) lhs <- est$lhs
  if(!('effect' %in% colnames(alpha))) {
    if(is.null(lhs))
        stop('Please specify lhs=[one of ',paste(est$lhs, collapse=','),']')
    effnam <- paste('effect',lhs,sep='.')
    if(!(effnam %in% colnames(alpha))) {
      stop("Can't find effect-column in alpha")
    }
  }
  resid <- est$residuals[,lhs]
  d1 <- alpha[alpha['fe']==nf1,effnam][f1]
  d2 <- alpha[alpha['fe']==nf2,effnam][f2]
  w <- if(is.null(est$weights)) 1.0 else est$weights
  if(is.null(est$weights)) {
    var1 <- var(d1)
    var2 <- var(d2)
    cov12 <- cov(d1,d2)
  } else {
    var1 <- wvar(d1,w^2)
    var2 <- wvar(d2,w^2)
    cov12 <- wcov(d1,d2,w^2)
  }
  delta1 <- varbias(corrfactors[[1]],est,tol,var1,maxsamples,resid=resid,weights=est$weights)
  if(is.null(est$weights))
      epsvar <- sum(resid^2)/est$df
  else
      epsvar <- sum(w^2*resid^2)/est$df/sum(w^2)

  if(nocovar) {
    f1 <- est$fe[[corrfactors[[1]]]]
    f2 <- est$fe[[corrfactors[[2]]]]
    N <- length(f1)
    delta2 <- delta1 - epsvar*(nlevels(f1)-nlevels(f2))/N
    delta12 <- epsvar*nlevels(f1)/N - delta1
  } else {
    delta2 <- varbias(corrfactors[[2]],est,tol,var2,maxsamples, resid=resid,weights=est$weights)
    eps <- -tol*sqrt((var1-delta1)*(var2-delta2))
    delta12 <- covbias(corrfactors,est,eps,maxsamples=maxsamples,resid=resid, weights=est$weights)
  }

  vartheta <- var1 - delta1
  varpsi <- var2 - delta2
  covtp <- cov12 - delta12

  c(corr=covtp/sqrt(vartheta*varpsi),
    v1=vartheta,
    v2=varpsi,
    cov=covtp,
    d1=delta1,
    d2=delta2,
    d12=delta12)

}

halftrace <- function(x, f, restf, MFX, tol, lmean, name, weights=NULL) {
  w <- weights
  if(is.null(w)) {ww <- ww2 <- 1.0} else {ww2 <- w^2; ww <- w}

  if(is.null(MFX)) {
    invfun <- function(v) {
      Crowsum(ww2*demeanlist(v[f,], restf, weights=w), f)
    }
  } else {
    invfun <- function(v) {
      Crowsum(ww2*demeanlist(demeanlist(ww*v[f,],MFX)/ww, restf, weights=w), f)
    }
  }

  DtM1x <- Crowsum(ww*demeanlist(x,lmean, weights=weights, scale=FALSE), f)
  # we use absolute tolerance, mctrace wil give us a trtol.
  # however, we square our result:  (Rx + err)^2 = Rx^2 + 2*err*Rx + err^2
  # so to get 2*err*Rx < tol, we must have err < tol/(2*Rx), i.e.
  # we should use a relative tolerance
#  tol1 <- -tol/sqrt(colSums(DtM1x^2))
  tol1 <- tol
#  message('cgsolve: tol1 ',tol1, ' tol ',tol)
  v <- cgsolve(invfun, DtM1x, eps=tol1,name=name)
 if(!is.null(MFX))
     Rx <- ww2*demeanlist(demeanlist(ww*v[f,],MFX)/ww, restf, weights=w)
 else
     Rx <- ww*demeanlist(v[f,],restf, weights=weights)

#  message("|v| = ",mean(sqrt(colSums(v^2))), '|Rx|=',mean(sqrt(colSums(Rx^2))))
  Rx
}

# In case we're weighted, the weighted variance of x is
# c  x' W M_W W x, not x' M_1 x
# W is the square root of the weights, c is sum(w)/(sum(w)^2 - sum(w^2))
# (where w is the weights (not the square root))
# that is, the bias is similar as before, but use M_W instead of M_1
# and Wx instead of x (i.e. WDtheta instead of Dtheta).
# remember that the weights in the felm-object are the square roots.
varbias <- function(index,est,tol=0.01,bvar, maxsamples=Inf,
                    robust=!is.null(est$clustervar), resid, weights=NULL,dfadj=1) {
  if(length(index) != 1) stop("index must have length 1")
  if(!is.null(est$stage1)) stop("Bias correction with IV not supported yet")
#  on.exit({rm(list=ls()); gc()})
  if(length(tol)==1) {
    tracetol <- tol
    cgtf <- 1
  } else {
    tracetol <- tol[[1]]
    cgtf <- tol[[2]]
  }
  X <- est$X
  f <- est$fe[[index]]
  restf <- est$fe[-index]
  name <- paste('var(',names(est$fe)[[index]],')', sep='')
  N <- length(f)
  nlev <- nlevels(f)
  w <- weights
  if(is.null(w)) {
    wc <- ww <- ww2 <- 1.0
  } else {
    w <- w/sqrt(sum(w^2))
    ww2 <- w^2; ww <- w; 
    wc <- 1/(1 - sum(ww2^2))
  }
  
# First, make a factor list for projecting out mean.
  fmean <- factor(rep(1,N))
  lmean <- list(fmean)

# project out X and F, i.e. everything but f1
# use M_{F,X} = M_F M_{M_F X}
# precompute and orthonormalize M_F X

  if(length(X) > 0) {
    # Find NA's in the coefficients, remove corresponding variables from X
    bad <- apply(est$coefficients,1,anyNA)
    if(any(bad)) X <- X[,!bad,drop=FALSE]
    MFX <- list(structure(fmean,x=orthonormalize(demeanlist(X,restf,
                                    weights=w,scale=c(TRUE,FALSE)))))
    invfun <- function(v) {
      Crowsum(ww*demeanlist(demeanlist(ww*v[f,],MFX), restf, weights=w, scale=FALSE), f)
    }
  } else {
    MFX <- NULL
    invfun <- function(v) {
      Crowsum(ww*demeanlist(v[f,], restf, weights=w,scale=c(TRUE,FALSE)), f)
    }
  }

  if(is.null(w))
      epsvar <- sum(resid^2)/est$df
  else
      epsvar <- sum(ww2*resid^2)*N/est$df

  vvfoo <- epsvar
  if(robust) {
    # residuals present, do cluster robust correction
    # create all the cluster interactions, so we don't have to do
    # it each time in the iteration
    docluster <- !is.null(est$clustervar)
    toladj <- sqrt(sum((ww*resid)^2))
    if(docluster) {
      d <- length(est$clustervar)
      cia <- list()
      for(i in 1:(2^d-1)) {
        # Find out which ones to interact
        iac <- as.logical(intToBits(i))[1:d]
        # interact the factors
        cia[[i]] <- factor(do.call(paste,c(est$clustervar[iac],sep='\004')))
      }
    }
    wwres <- ww*resid
    trfun <- function(x,trtol) {
      # return crude estimate of the trace
      if(trtol == 0) return(abs(nlev))
#      on.exit({rm(list=ls()); gc()})
      # since we square our result, we should use sqrt(trtol)/2 as a tolerance
      # trtol is the absolute tolerance, but halftrace is squared, we should
      # really use a relative tolerance, we use the square root of the tracetol
#      message('trtol=',trtol, ' N ',N,' toladj ',toladj)
      Rx <- halftrace(x, f, restf, MFX, -sqrt(trtol/toladj), lmean, name, weights=w)
#      Rx <- halftrace(x, f, restf, MFX, sqrt(trtol)/2, lmean, name, weights=w)
#      Rx <- halftrace(x, f, restf, MFX, trtol, lmean, name, weights=w)

      # now apply cluster stuff
      # first, scale with (weighted) residuals
#
      .Call(C_scalecols, Rx, wwres)
      if(!docluster) {
        # It's heteroscedastic
        return(colSums(Rx * Rx))
      } else {
        # it's one or more clusters, do the Cameron et al detour
        result <- vector('numeric',ncol(Rx))
        for(i in 1:(2^d-1)) {
          ia <- cia[[i]]
          b <- Crowsum(Rx,ia)
          # odd number is positive, even is negative
          sgn <- 2*(sum(as.logical(intToBits(i))[1:d]) %% 2) - 1
          adj <- sgn*dfadj*nlevels(ia)/(nlevels(ia)-1)
          result <- result + adj* colSums(b * b)
          rm(b)
#          result <- result + vvfoo*adj* colSums(Rx * Rx)
        }
        return(result)
      }
    }
    epsvar <- 1  # now incorporated in the robust variance matrix, so don't scale
  } else {

    trfun <- function(x,trtol) {
      # return crude estimate of the trace
      if(trtol == 0) return(abs(nlev))
#      on.exit({rm(list=ls()); gc()})
      DtM1x <- Crowsum(ww*demeanlist(unnamed(x),lmean,weights=w,scale=FALSE), f)
      # we use absolute tolerance, mctrace wil give us a trtol.
      # we divide by the L2-norm of DtM1x, since we take the
      # inner product with this afterwards

      tol1 <- -trtol/cgtf/sqrt(colSums(DtM1x^2))/2
      v <- cgsolve(invfun, DtM1x, eps=tol1,name=name)
      colSums(DtM1x * v)
    }
  }
  attr(trfun,'IP') <- TRUE
  epsvar <- epsvar * wc

  # If we have weights, we should not divide the trace by N

  # epsvar*trace/N is the bias estimate
  # We want precision in the final estimate to be, say, 1%
  # Since we subtract the bias from the biased estimate, the precision
  # in the trace computation depends on the current value of the trace
  # i.e. absolute precision should be 0.01*(bvar*N - epsvar*tr)/epsvar
  # where bvar is the biased variance

  epsfun <- function(tr) {
    aa <- N*bvar - epsvar*tr
    if(aa < 0) {
      -abs(tracetol)*0.5*N*bvar/epsvar
    } else {
      -abs(tracetol)*aa/epsvar
    }
  }
  # the tolerance before mctrace has got a clue about where we are,
  # is a problem. If the bias is very large compared to the variance, we will
  # be in trouble. 
  res <- epsvar*mctrace(trfun,N=N,tol=epsfun, trname=name,
                 maxsamples=maxsamples)/N
}


# if positive tolerance, the tolerance is relative to the bias corrected
# covariance. In this case, the biased covariance (bcov) and residual
# variance (epsvar) must be specified. A negative tolerance is an
# absolute tolerance
covbias <- function(index,est,tol=0.01, maxsamples=Inf, resid, weights=NULL,
                    robust=!is.null(est$clustervar)) {
  if(length(index) != 2) stop("index must have length 2")
  if(!is.null(est$stage1)) stop("Bias correction with IV not supported yet")
  if(length(est$fe) < 2) stop("fe must have length >= 2")
#  on.exit({rm(list=ls()); gc()})
  w <- weights
  if(is.null(w)) {
    wc <- ww2 <- ww <- 1
  } else {
    w <- w/sqrt(sum(w^2))
    ww <- w; ww2 <- w^2
    wc <- sum(w^2)/(sum(w^2)^2 - sum(w^4))
  }
  X <- est$X
  f1 <- est$fe[[index[[1]]]]
  f2 <- est$fe[[index[[2]]]]
  nlev1 <- nlevels(f1)
  nlev2 <- nlevels(f2)
  N <- length(f1)
  name <- paste('cov(',paste(names(est$fe)[index],collapse=','),')',sep='')
  name1 <- paste(name,names(est$fe)[index[[1]]],sep='.')
  name2 <- paste(name,names(est$fe)[index[[2]]], sep='.')
  no2list <- est$fe[-index[[2]]]
  no1list <- est$fe[-index[[1]]]

  restf <- est$fe[-index]
  fmean <- factor(rep(1,N))
  lmean <- list(fmean)
  if(is.null(w))
      epsvar <- sum(resid^2)/est$df
  else
      epsvar <- sum(ww2*resid^2)*N/est$df

  if(length(X) > 0) {
    bad <- apply(est$coefficients,1,anyNA)
    if(any(bad)) X <- X[,!bad,drop=FALSE]
  }
  MDX <- MFX <- NULL

  if(length(X) > 0) {
    MDX <- list(structure(fmean,x=orthonormalize(demeanlist(X,
                                    no2list,weights=w, scale=c(TRUE,FALSE)))))
    MFX <- list(structure(fmean,x=orthonormalize(demeanlist(X,
                                    no1list,weights=w, scale=c(TRUE,FALSE)))))

    invfun <- function(v) {
      Crowsum(ww*demeanlist(demeanlist(ww*v[f2,],MDX), no2list,
                           weights=w, scale=FALSE), f2)
    }
    if(length(restf) > 0) {
      MX <- list(structure(fmean,x=orthonormalize(demeanlist(X,
                                     restf,weights=w, scale=c(TRUE,FALSE)))))
      MXfun <- function(v) ww*demeanlist(demeanlist(ww*v[f1,], MX), restf, weights=w,
                                         scale=FALSE)
    } else {
      MX <- list(structure(fmean, x=orthonormalize(ww*X)))
      MXfun <- function(v) ww*demeanlist(ww*v[f1,], MX)
    }
  } else {
    invfun <- function(v) {
      Crowsum(ww2*demeanlist(v[f2,], no2list, weights=w), f2)
    }
    MXfun <- function(v) ww2*demeanlist(v[f1,], restf, weights=w)
  }

  invfunX <- function(v) {
    Crowsum(MXfun(v), f1)
  }
  
  if(robust) {
    # residuals present, do cluster robust correction
    # create all the cluster interactions, so we don't have to do
    # it each time in the iteration
    tolmod <- 1
    docluster <- !is.null(est$clustervar)
    if(docluster) {
      tolmod <- 0
      d <- length(est$clustervar)
      cia <- list()
      for(i in 1:(2^d-1)) {
        # Find out which ones to interact
        iac <- as.logical(intToBits(i))[1:d]
        # interact the factors
        cia[[i]] <- factor(do.call(paste,c(est$clustervar[iac],sep='\004')))
        tolmod <- tolmod + nlevels(cia[[i]])
      }
    }
    toladj <- sqrt(sum((ww*resid)^2))
    wwres <- ww*resid
    trfun <- function(x,trtol) {
      # return crude estimate of the trace
      if(trtol == 0) return(abs(nlev1-nlev2))
#      on.exit({rm(list=ls()); gc()})
      # since we square our result, we should use sqrt(trtol)/2 as a tolerance
    # trtol is the absolute tolerance, but halftrace is squared, we should
      # really use a relative tolerance, we use the square root of the tracetol
#      message('cov trtol=',trtol, ' tolmod ',tolmod, ' N ',N)
      Lx <- halftrace(x, f1, no1list, MFX, -sqrt(trtol/toladj), lmean, name1, weights=w)
      Rx <- halftrace(x, f2, no2list, MDX, -sqrt(trtol/toladj), lmean, name2, weights=w)
#      Rx <- halftrace(x, f, restf, MFX, sqrt(trtol)/2, lmean, name, weights=w)
#      Rx <- halftrace(x, f, restf, MFX, trtol, lmean, name, weights=w)

      # now apply cluster stuff
      # first, scale with (weighted) residuals
      .Call(C_scalecols, Rx, wwres)
      .Call(C_scalecols, Lx, wwres)
      if(!docluster) {
        # It's heteroscedastic
        return(colSums(Lx * Rx))
      } else {
        # it's one or more clusters, do the Cameron et al detour
        result <- vector('numeric',ncol(Rx))
        for(i in 1:(2^d-1)) {
          Lb <- Crowsum(Lx, cia[[i]])
          Rb <- Crowsum(Rx, cia[[i]])
          # odd number is positive, even is negative
          sgn <- 2*(sum(as.logical(intToBits(i))[1:d]) %% 2) - 1
          result <- result + sgn * colSums(Lb * Rb)
        }
        return(result)
      }
    }
    epsvar <- 1  # now incorporated in the robust variance matrix, so don't scale
  } else {
    trfun <- function(x,trtol) {
      # return crude estimate of the trace
      if(trtol == 0) return(-abs(nlev1-nlev2))
#      on.exit({rm(list=ls()); gc()})
      M1x <- ww*demeanlist(unnamed(x),lmean,weights=w,scale=FALSE)
      DtM1x <- Crowsum(M1x,f1)
      FtM1x <- Crowsum(M1x,f2)
      d1 <- colSums(DtM1x^2)
      d2 <- colSums(FtM1x^2)
      v <- cgsolve(invfunX, DtM1x, eps=-trtol/sqrt(d1+d2)/3, name=name)
      MXv <- Crowsum(MXfun(v), f2)
      sol <- cgsolve(invfun, FtM1x, eps=-trtol/sqrt(colSums(MXv^2))/3, name=name)
      -colSums(sol* MXv)
    }
  }
  # our function does the inner product, not just matrix application. Signal to mctrace.
  attr(trfun,'IP') <- TRUE
  epsvar <- epsvar * wc

  # absolute precision, scale by N and epsvar since it's trace level precision

  eps <- -abs(tol)*N/epsvar

  epsvar*mctrace(trfun, N=N, tol=eps, trname=name,
                 maxsamples=maxsamples)/N
}

# compute variance of the biased variance estimate
# using the formula for the variance of a quadratic form with normal distribution
# var(x^t A x) = 2 tr(AVAV) + 4mu^t*AVA*mu
# where V is the variance matrix of x, assumed to be sigma^2 I, and mu is the
# expectation of x (i.e. Dtheta).
varvar <- function(index, fe, X, pointest, resvar, tol=0.01,
                   biascorrect=FALSE, weights=NULL) {
  w <- weights
#  w <- NULL
  if(is.null(w)) {
    wc <- ww <- ww2 <- 1.0
  } else {
    w <- w/sqrt(sum(w^2))
    ww2 <- w^2; ww <- w
    wc <- sum(w^2)/(sum(w^2)^2 - sum(w^4))
  }
  if(!is.null(w) && biascorrect) warning('bias corrected varvars with weights not tested')
  f <- fe[[index]]
  N <- length(f)
  lmean <- list(factor(rep(1,N)))
  name <- paste('varvar(',names(fe)[[index]],')', sep='')
  if(length(X)==0) {
      MFX <- fe[-index]
      invfun <- function(x) {
        Crowsum(ww2*demeanlist(x[f,], MFX, weights=w),f)
      }
  } else {
#    M_{F,X} = M_F M_{M_F X}
    restf <- fe[-index]
    MFX <- list(structure(factor(rep(1,N)),
                          x=orthonormalize(ww*demeanlist(X,restf,weights=w))))
    invfun <- function(x) {
      Crowsum(ww2*demeanlist(demeanlist(ww*x[f,],MFX)/ww, restf, weights=w), f)
    }
  }


  Dtheta <- pointest[f]
  DtM1D <- Crowsum(ww2*demeanlist(Dtheta,lmean,weights=w), f)
  v <- cgsolve(invfun, DtM1D, eps=tol/4/resvar/sqrt(sum(DtM1D^2)), name=name)
  meanpart <- 4*resvar * sum(DtM1D * v)
  if(!biascorrect) return(meanpart/N^2)
  # the mean part is biased upwards. We should correct it.
  # it turns out that we can do this by changing the sign of the
  # trace term, the bias is the same expression as the trace part
  #  message('mean part=',meanpart/N^2)
  mytol <- meanpart/10
  trfun <- function(x,trtol) {
    v <- ww*demeanlist(cgsolve(invfun, Crowsum(ww2*demeanlist(x,lmean,weights=w),f),
                            eps=-mytol^2/resvar^2/2,name=name)[f,],lmean,weights=w)
    colSums(v * x)
  }
  attr(trfun,'IP') = TRUE
  trpart <- 2*resvar^2 * mctrace(trfun, N=length(f), trname=name, tol=-mytol/resvar/2)
  if(!is.null(w)) trpart <- trpart/N
#  message('mean part=', meanpart, ' trpart=',trpart)
  (meanpart-trpart)/N^2
}




#' Compute the variance of the fixed effect variance estimate
#'
#' @details
#' With a model like \eqn{y = X\beta + D\theta + F\psi + \epsilon}, where \eqn{D} and
#' \eqn{F} are matrices with dummy encoded factors, one application of \pkg{lfe} is
#' to study the variances \eqn{var(D\theta)}, \eqn{var(F\psi)} and covariances
#' \eqn{cov(D\theta, F\psi)}. The function \code{\link{fevcov}} computes bias corrected
#' variances and covariances.  However, these variance estimates are still
#' random variables for which \code{\link{fevcov}} only estimate the
#' expectation. The function \code{varvars} estimates the variance of these
#' estimates.
#' 
#' This function returns valid results only for normally distributed residuals.
#' Note that the estimates for the fixed effect variances from
#' \code{\link{fevcov}} are not normally distributed, but a sum of chi-square
#' distributions which depends on the eigenvalues of certain large matrices. We
#' do not compute that distribution. The variances returned by \code{varvars}
#' can therefore \emph{not} be used directly to estimate confidence intervals,
#' other than through coarse methods like the Chebyshev inequality. These
#' estimates only serve as a rough guideline as to how wrong the variance
#' estimates from \code{\link{fevcov}} might be.
#' 
#' Like the fixed effect variances themselves, their variances are also biased
#' upwards.  Correcting this bias can be costly, and is therefore by default
#' switched off.
#' 
#' The variances tend to zero with increasing number of observations. Thus, for
#' large datasets they will be quite small.
#' 
#' @param est an object of class '"felm"', the result of a call to
#' \code{\link{felm}(keepX=TRUE)}.
#' @param alpha a data frame, the result of a call to \code{\link{getfe}}.
#' @param tol numeric. The absolute tolerance for the bias-corrected
#' correlation.
#' @param biascorrect logical. Should the estimates be bias corrected?
#' @param lhs character. Name of left hand side if multiple left hand sides.
#' @return \code{varvars} returns a vector with a variance estimate for each
#' fixed effect variance.  I.e. for the diagonal returned by
#' \code{\link{fevcov}}.
#' @note The \code{tol} argument specifies the tolerance as in
#' \code{\link{fevcov}}.  Note that if \code{est} is the result of a call to
#' \code{\link{felm}} with \code{keepX=FALSE} (the default), the variances will
#' be estimated as if the covariates X are independent of the factors.  There
#' is currently no function available for estimating the variance of the
#' covariance estimates from \code{\link{fevcov}}.
#' 
#' The cited paper does not contain the expressions for the variances computed
#' by \code{varvars} (there's a 10 page limit in that journal), though they can
#' be derived in the same fashion as in the paper, with the formula for the
#' variance of a quadratic form.
#' @seealso \code{\link{bccorr}} \code{\link{fevcov}}
#' @references Gaure, S. (2014), \cite{Correlation bias correction in two-way
#' fixed-effects linear regression}, Stat 3(1):379-390, 2014.
#' @examples
#' 
#' x <- rnorm(500)
#' x2 <- rnorm(length(x))
#' 
#' ## create individual and firm
#' id <- factor(sample(40,length(x),replace=TRUE))
#' firm <- factor(sample(30,length(x),replace=TRUE,prob=c(2,rep(1,29))))
#' foo <- factor(sample(20,length(x),replace=TRUE))
#' ## effects
#' id.eff <- rnorm(nlevels(id))
#' firm.eff <- rnorm(nlevels(firm))
#' foo.eff <- rnorm(nlevels(foo))
#' ## left hand side
#' id.m <- id.eff[id]
#' firm.m <- 2*firm.eff[firm]
#' foo.m <- 3*foo.eff[foo]
#' y <- x + 0.25*x2 + id.m + firm.m + foo.m + rnorm(length(x))
#' 
#' # make a data frame
#' fr <- data.frame(y,x,x2,id,firm,foo)
#' ## estimate and print result
#' est <- felm(y ~ x+x2|id+firm+foo, data=fr, keepX=TRUE)
#' alpha <- getfe(est)
#' # estimate the covariance matrix of the fixed effects
#' fevcov(est, alpha)
#' # estimate variances of the diagonal
#' varvars(est, alpha)
#' 
#' @export varvars
varvars <- function(est, alpha=getfe(est), tol=0.01, biascorrect=FALSE, lhs=NULL) {
  if(nlevels(est$cfactor) > 1) stop('Data should have just a single connected component')

  fe <- est$fe
  e <- length(fe)
  if(length(tol) == 1) tol <- rep(tol,e)

  if(length(est$lhs) > 1 && is.null(lhs))
      stop('Please specify lhs=[one of ',paste(est$lhs, collapse=','),']')      
  if(length(est$lhs) == 1) lhs <- est$lhs
  effnam <- 'effect'
  if(! ('effect' %in% colnames(alpha))) {
    effnam <- paste('effect',lhs,sep='.')
    if(!(effnam %in% colnames(alpha))) {
      stop("Can't find effect-column in alpha")
    }
  }

  effs <- lapply(names(fe), function(nm) alpha[alpha[,'fe']==nm, effnam])
  w2 <- if(is.null(est$weights)) 1 else est$weights^2
  resvar <- sum(w2*est$residuals[,lhs]^2)*length(w2)/est$df

  sapply(1:e, function(index) {
    varvar(index, fe, est$X, effs[[index]], resvar, tol[index], biascorrect, est$weights)
  })
}



# a function for computing the covariance matrix between
# all the fixed effects






#' Compute limited mobility bias corrected covariance matrix between fixed
#' effects
#' 
#' With a model like \eqn{y = X\beta + D\theta + F\psi + \epsilon}, where
#' \eqn{D} and \eqn{F} are matrices with dummy encoded factors, one application
#' of \pkg{lfe} is to study the variances \eqn{var(D\theta)}, \eqn{var(F\psi)}
#' and covariances \eqn{cov(D\theta, F\psi)}. However, if we use estimates for
#' \eqn{\theta} and \eqn{\psi}, the resulting variances are biased. The
#' function \code{fevcov} computes a bias corrected covariance matrix as
#' described in \cite{Gaure (2014)}.
#' 
#' The \code{tol} argument specifies the tolerance. The tolerance is relative
#' for the variances, i.e. the diagonal of the output.  For the covariances,
#' the tolerance is relative to the square root of the product of the
#' variances, i.e. an absolute tolerance for the correlation.  If a numeric of
#' length 1, \code{tol} specifies the same tolerance for all
#' variances/covariances.  If it is of length 2, \code{tol[1]} specifies the
#' variance tolerance, and \code{tol[2]} the covariance tolerance.  \code{tol}
#' can also be a square matrix of size \code{length(est$fe)}, in which case the
#' tolerance for each variance and covariance is specified individually.
#' 
#' The function performs no checks for estimability. If the fixed effects are
#' not estimable, the result of a call to \code{fevcov} is not useable.
#' Moreover, there should be just a single connected component among the fixed
#' effects.
#' 
#' \code{alpha} must contain a full set of coefficients, and contain columns
#' \code{'fe'} and \code{'effect'} like the default estimable functions from
#' \code{\link{efactory}}.
#' 
#' In the case that the \code{\link{felm}}-estimation has weights, it is the
#' weighted variances and covariance which are bias corrected.
#' 
#' @param est an object of class '"felm"', the result of a call to
#' \code{\link{felm}(keepX=TRUE)}.
#' @param alpha a data frame, the result of a call to \code{\link{getfe}}.
#' @param tol numeric. The absolute tolerance for the bias-corrected
#' correlation.
#' @param robust logical. Should robust (heteroskedastic or cluster) residuals
#' be used, rather than i.i.d.
#' @param maxsamples integer. Maximum number of samples for expectation
#' estimates.
#' @param lhs character. Name of left hand side if multiple left hand sides.
#' @return \code{fevcov} returns a square matrix with the bias corrected
#' covariances. An attribute \code{'bias'} contains the biases.  The bias
#' corrections have been subtracted from the bias estimates.  I.e. vc = vc' -
#' b, where vc' is the biased variance and b is the bias.
#' @note Bias correction for IV-estimates are not supported as of now.
#' 
#' Note that if \code{est} is the result of a call to \code{\link{felm}} with
#' \code{keepX=FALSE} (the default), the biases will be computed as if the
#' covariates X are independent of the factors. This will be faster (typically
#' by a factor of approx. 4), and possibly wronger.  Note also that the
#' computations performed by this function are non-trivial, they may take quite
#' some time.  It would be wise to start out with quite liberal tolerances,
#' e.g. \cite{tol=0.1}, to get an idea of the time requirements.
#' 
#' If there are only two fixed effects, \code{fevcov} returns the same
#' information as \code{\link{bccorr}}, though in a slightly different format.
#' @seealso \code{\link{varvars}} \code{\link{bccorr}}
#' @references Gaure, S. (2014), \cite{Correlation bias correction in two-way
#' fixed-effects linear regression}, Stat 3(1):379-390, 2014.
#' \url{http://dx.doi.org/10.1002/sta4.68}
#' @examples
#' 
#' x <- rnorm(5000)
#' x2 <- rnorm(length(x))
#' 
#' ## create individual and firm
#' id <- factor(sample(40,length(x),replace=TRUE))
#' firm <- factor(sample(30,length(x),replace=TRUE,prob=c(2,rep(1,29))))
#' foo <- factor(sample(20,length(x),replace=TRUE))
#' ## effects
#' id.eff <- rnorm(nlevels(id))
#' firm.eff <- runif(nlevels(firm))
#' foo.eff <- rchisq(nlevels(foo),df=1)
#' ## left hand side
#' id.m <- id.eff[id]
#' firm.m <- firm.eff[firm]
#' foo.m <- foo.eff[foo]
#' # normalize them
#' id.m <- id.m/sd(id.m)
#' firm.m <- firm.m/sd(firm.m)
#' foo.m <- foo.m/sd(foo.m)
#' y <- x + 0.25*x2 + id.m + firm.m + foo.m + rnorm(length(x),sd=2)
#' z <- x + 0.5*x2 + 0.7*id.m + 0.5*firm.m + 0.3*foo.m + rnorm(length(x),sd=2)
#' # make a data frame
#' fr <- data.frame(y,z,x,x2,id,firm,foo)
#' ## estimate and print result
#' est <- felm(y|z ~ x+x2|id+firm+foo, data=fr, keepX=TRUE)
#' # find bias corrections, there's little bias in this example
#' print(yv <- fevcov(est, lhs='y'))
#' ## Here's how to compute the unbiased correlation matrix:
#' cm <- cov2cor(yv)
#' structure(cm,bias=NULL)
#' 
#' @export fevcov
fevcov <- function(est, alpha=getfe(est), tol=0.01, robust=!is.null(est$clustervar),
                   maxsamples=Inf, lhs=NULL) {
  if(nlevels(est$cfactor) > 1) stop('Data should have just a single connected component')
  iaf <- sapply(est$fe, function(f) !is.null(attr(f,'x',exact=TRUE)))
  if(any(iaf))
      stop("Bias correction for interacted factor ",paste(names(est$fe)[iaf],collapse=', '), " not supported")
  if(length(est$lhs) > 1 && is.null(lhs))
      stop('Please specify lhs=[one of ',paste(est$lhs, collapse=','),']')      
  if(length(est$lhs) == 1) lhs <- est$lhs
  if(is.null(lhs)) lhs <- 1
  effnam <- 'effect'
  if(! ('effect' %in% colnames(alpha))) {
    effnam <- paste('effect',lhs,sep='.')
    if(!(effnam %in% colnames(alpha))) {
      stop("Can't find effect-column in alpha")
    }
  }
  if(is.na(match('fe', colnames(alpha))))
      stop("alpha should contain columns 'fe' and 'effect'")

#  if(!is.null(est$weights)) warning("Bias correction with weights not yet fully correct")

  if(length(tol) == 1) tol <- c(tol,-abs(tol))
  if(length(tol) != 2 && !is.matrix(tol))
      stop('tol must be either a matrix or of length 1 or 2')
  K <- length(est$fe)
  if(is.matrix(tol) && (ncol(tol) != K || nrow(tol) != K))
      stop('tol matrix must be square, with size the number of fixed effects: ',K)
  if(!is.matrix(tol)) {
    tmp <- tol
    tol <- matrix(0,K,K)
    diag(tol) <- tmp[1]
    tol[col(tol) != row(tol)] <- tmp[2]
  }

  # compute the biased variances
  fe <- est$fe
  effs <- lapply(names(fe), function(nm) alpha[alpha[,'fe']==nm, effnam][fe[[nm]]])
  names(effs) <- names(fe)
  bvcv <- matrix(0,K,K)
  colnames(bvcv) <- rownames(bvcv) <- names(fe)
  diag(bvcv) <- sapply(effs, var)
  for(i in 1:K) {
    for(j in i:K) {
      if(is.null(est$weights)) {
        bvcv[i,j] <- bvcv[j,i] <- cov(effs[[i]],effs[[j]])
      } else {
        bvcv[i,j] <- bvcv[j,i] <- wcov(effs[[i]], effs[[j]], est$weights^2)
      }
    }
  }
  # compute the variances
  bias <- matrix(0,K,K)
  colnames(bias) <- rownames(bias) <- names(fe)
  resid <- est$residuals[,lhs]
  cgtf <- rep(1,K)
  diag(bias) <- sapply(1:K, function(i) varbias(i, est, tol[i,i], bvcv[i,i],
                                                maxsamples, resid=resid,weights=est$weights,
                                                dfadj=(est$N-1)/est$df))
  
  wtf <- diag(bias) > diag(bvcv)
  # update off-diagonal tolerances by the variances
  offdiag <- col(tol) != row(tol)
  if(any(wtf)) {
    message('Some variance biases are larger than the variances. Setting them equal:')
    print(cbind(variance=diag(bvcv),bias=diag(bias)))
    diag(bias)[wtf] <- 0.999*diag(bvcv)[wtf]
    tol[offdiag] <- -(abs(tol)*sqrt(abs(tcrossprod(diag(bvcv)-0.9*diag(bias)))))[offdiag]
  } else {
    tol[offdiag] <- -(abs(tol)*sqrt(abs(tcrossprod(diag(bvcv)-diag(bias)))))[offdiag]
  }
  # compute the covariances
  if(K > 1) {
    for(i in 1:(K-1)) {
      for(j in (i+1):K)
          bias[i,j] <- bias[j,i] <- covbias(c(i,j),est,tol[i,j],maxsamples,
                                            resid=resid, weights=est$weights)
    }  
  }
  structure(bvcv-bias, bias=bias)
}
