#' Bootstrap standard errors for the group fixed effects
#' 
#' Bootstrap standard errors for the group fixed effects which were swept out
#' during an estimation with \code{\link{felm}}.
#' 
#' The bootstrapping is done in parallel if \code{threads > 1}.
#' \code{\link{btrap}} is run automatically from \code{\link{getfe}} if
#' \code{se=TRUE} is specified.  To save some overhead, the individual
#' iterations are grouped together, the memory available for this grouping is
#' fetched with \code{getOption('lfe.bootmem')}, which is initialized upon
#' loading of \pkg{lfe} to \code{options(lfe.bootmem=500)} (MB).
#' 
#' If \code{robust=TRUE}, heteroskedastic robust standard errors are estimated.
#' If \code{robust=FALSE} and \code{cluster=TRUE}, clustered standard errors
#' with the cluster specified to \code{felm()} are estimated. If \code{cluster}
#' is a factor, it is used for the cluster definition.  \code{cluster may} also
#' be a list of factors.
#' 
#' @param alpha data frame returned from \code{\link{getfe}}
#' @param obj object of class \code{"felm"}, usually, a result of a call to
#' \code{\link{felm}}
#' @param N integer.  The number of bootstrap iterations
#' @param ef function.  An estimable function such as in \code{\link{getfe}}.
#' The default is to use the one used on \code{alpha}
#' @param eps double. Tolerance for centering, as in getfe
#' @param threads integer.  The number of threads to use
#' @param robust logical. Should heteroskedastic standard errors be estimated?
#' @param cluster logical or factor. Estimate clustered standard errors.
#' @param lhs character vector. Specify which left hand side if \code{obj} has
#' multiple lhs.
#' @return A data-frame of the same size as alpha is returned, with standard
#' errors filled in.
#' @examples
#' 
#' oldopts <- options(lfe.threads=2)
#' ## create covariates
#' x <- rnorm(3000)
#' x2 <- rnorm(length(x))
#' 
#' ## create individual and firm
#' id <- factor(sample(700,length(x),replace=TRUE))
#' firm <- factor(sample(300,length(x),replace=TRUE))
#' 
#' ## effects
#' id.eff <- rlnorm(nlevels(id))
#' firm.eff <- rexp(nlevels(firm))
#' 
#' ## left hand side
#' y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))
#' 
#' ## estimate and print result
#' est <- felm(y ~ x+x2 | id + firm)
#' summary(est)
#' ## extract the group effects
#' alpha <- getfe(est)
#' head(alpha)
#' ## bootstrap standard errors
#' head(btrap(alpha,est))
#' 
#' ## bootstrap some differences
#' ef <- function(v,addnames) {
#'   w <- c(v[2]-v[1],v[3]-v[2],v[3]-v[1])
#'   if(addnames) {
#'      names(w) <-c('id2-id1','id3-id2','id3-id1')
#'      attr(w,'extra') <- list(note=c('line1','line2','line3'))
#'   }
#'   w
#' }
#' # check that it's estimable
#' is.estimable(ef,est$fe)
#' 
#' head(btrap(alpha,est,ef=ef))
#' options(oldopts)
#' 
#' @export btrap
btrap <- function(alpha,obj,N=100,ef=NULL,eps=getOption('lfe.eps'),
                  threads=getOption('lfe.threads'), robust=FALSE,
                  cluster=NULL, lhs=NULL) {
  # bootstrap the stuff. The name 'btrap' is chosen to instill a feeling of being trapped
  # (in some long running stuff which will never complete)
  # bootstrapping is really to draw residuals over again, i.e. to change
  # the outcome.  Predictions of the estimated system are adjusted by
  # drawing from the residuals. We need a new r.h.s to replace
  # the current (I-P)(Y-Xbeta), i.e. (I-P)(Y-Xbeta+delta) where delta is resampled from
  # the residuals PY-PXbeta.  Then we adjust for degrees of freedom, which is component specific.
    
  if(is.logical(cluster)) {
      if(cluster) cluster <- obj$clustervar else cluster <- NULL
  } else if(!is.null(cluster)) {
      if(is.list(cluster))
          cluster <- lapply(cluster,factor)
      else
          cluster <- list(factor(cluster))
  }

  if(is.null(ef))  {
    # use the one used with alpha
    ef <- attr(alpha,'ef')
  } else {
    if(is.character(ef)) ef <- efactory(obj,opt=ef)
    # redo the point estimates
    v <- ef(alpha[,'effect'],TRUE)
    alpha <- data.frame(effect=v)
    rownames(alpha) <- names(v)
    if(!is.null(attr(v,'extra'))) alpha <- cbind(alpha,attr(v,'extra'))
  }

  if(is.null(lhs)) {
    R <- obj$r.residuals-obj$residuals
    smpdraw <- as.vector(obj$residuals)
  } else {
    R <- obj$r.residuals[,lhs]-obj$residuals[,lhs]
    smpdraw <- as.vector(obj$residuals[,lhs])
  }

  w <- obj$weights
  # if there are weights, smpdraw should be weighted
  if(!is.null(w)) smpdraw <- smpdraw * w

  # Now, we want to do everything in parallel, so we should allocate up a set
  # of vectors, but we don't want to blow the memory.  Stick to allocating two
  # vectors per thread.  The threaded stuff can't be interrupted, so this is
  # an opportunity to control-c too.
  # hmm, up to 500 MB of vectors, we say, but no less than two per thread
  # (one per thread is bad for balance, if time to completion varies)
  # divide by two because we use a copy in the demeanlist step.
  maxB <- getOption('lfe.bootmem')*1e6/2
  vpt <- max(2,as.integer(min(maxB/(length(R)*8),N)/threads))
  vpb <- vpt*threads
  blks <- as.integer(ceiling(N / vpb))
  newN <- blks*vpb
  vsum <- 0
  vsq <- 0
  start <- last <- as.integer(Sys.time())
  gc()

  if(is.null(lhs))
      predy <- obj$fitted.values
  else
      predy <- obj$fitted.values[,lhs]

  if(!is.null(cluster)) {
  # now, what about multiway clustering?
    # try the Cameron-Gelbach-Miller stuff.
    # make the interacted factors, with sign
    iac <- list()
    d <- length(cluster)
    for(i in 1:(2^d-1)) {
      # Find out which ones to interact
      iab <- as.logical(intToBits(i))[1:d]
      # odd number is positive, even is negative
      sgn <- 2*(sum(iab) %% 2) - 1
      iac[[i]] <- list(sgn=sgn/choose(d,sum(iab)),ib=iab)
    }
  }

  X <- model.matrix(obj,centred=NA)
  cX <- attr(X,'cX')

  for(i in 1:blks) {
    if(robust) {
      # robust residuals, variance is each squared residual
      #      rsamp <- rnorm(vpb*length(smpdraw))*abs(smpdraw)
      rsamp <- lapply(1:vpb,function(i) rnorm(length(smpdraw))*abs(smpdraw))
    } else if(!is.null(cluster)) {
      # Wild bootstrap with Rademacher distribution
      rsamp <- lapply(1:vpb, function (i) {
        if(length(cluster) == 1) {
          smpdraw*sample(c(-1,1),nlevels(cluster[[1]]), replace=TRUE)[cluster[[1]]]
        } else {
          # draw a Rademacher dist for each single cluster
          rad <- lapply(cluster,function(f) sample(c(-1,1),nlevels(f),replace=TRUE)[f])
          Reduce('+',lapply(iac, function(ia) ia$sgn*smpdraw * Reduce('*',rad[ia$ib])))
        }
      })
    } else {
      # IID residuals
      rsamp <- lapply(1:vpb, function(i) sample(smpdraw, replace=TRUE))
    }
    if(!is.null(w)) rsamp <- lapply(rsamp, function(x) x/w)
    if(length(cX) > 0) {
      newr <- lapply(rsamp, function(rs) {
        newy <- predy + rs
        newbeta <- obj$inv %*% crossprod(cX, newy)
        newbeta[is.na(newbeta)] <- 0
        as.vector(newy - X %*% newbeta)
      })
    } else {
      newr <- lapply(rsamp, function(rs) as.vector(predy) + rs)
    }
    rm(rsamp)
    v <- kaczmarz(obj$fe, demeanlist(unnamed(newr), obj$fe, eps=eps, threads=threads,
                                     means=TRUE, weights=w),
                  eps=eps, threads=threads)
    rm(newr)
    efv <- lapply(v,ef,addnames=FALSE)
    vsum <- vsum + Reduce('+',efv)
    vsq <- vsq + Reduce('+',Map(function(i) i^2, efv))
    now <- as.integer(Sys.time())
    if(now-last > 300) {
      cat('...finished',i*vpb,'of',newN,'vectors in',now-start,'seconds\n')
      last <- now
    }
  }
  if(robust)
      sename <- 'robustse'
  else if(!is.null(cluster))
      sename <- 'clusterse'
  else
      sename <- 'se'
  fSEname <- SEname <- 'se'
  fsename <- sename
  if(!is.null(lhs)) {fsename <- paste(sename,lhs,sep='.'); fSEname <- paste(SEname,lhs,sep='.');}
  alpha[,fsename] <- sqrt(vsq/newN - (vsum/newN)**2)/(1-0.75/newN-7/32/newN**2-9/128/newN**3)
  if(sename != 'se') alpha[,fSEname] <- alpha[,fsename]
  return(structure(alpha,sename=sename))
}
