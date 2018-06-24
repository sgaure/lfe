# $Id: felm.old.R 1655 2015-03-18 18:51:06Z sgaure $
felm.old <- function(formula,fl,data) {
  mf <- match.call(expand.dots = FALSE)

  if(missing(fl)) {
    # we should rather parse the formula tree
    # find the terms involving G
    trm <- terms(formula,special='G')
    feidx <- attr(trm,'specials')$G-1
    festr <- paste(labels(trm)[feidx],collapse='+')

    if(festr == '') stop('No factors specified')
    # remove the G-terms from formula
    formula <- update(formula,paste('. ~ . -(',festr,')'))
    mf[['formula']] <- formula

    # then make a list of them, and find their names
    felist <- parse(text=paste('list(',gsub('+',',',festr,fixed=TRUE),')',sep=''))
    nm <- eval(felist,list(G=function(t) as.character(substitute(t))))
    # collapse them in case there's an interaction with a funny name
    nm <- lapply(nm,paste,collapse='.')
    
    # replace G with as.factor, eval with this, and the parent frame, or with data
    # allow interaction factors with '*'
    iact <- function(a,b) interaction(a,b,drop=TRUE)
    if(missing(data)) 
      fl <- eval(felist,list(G=as.factor,'*'=iact))
    else {
      G <- as.factor
      fl <- local({'*'<-iact;eval(felist,data,environment())})
    }
    gc()
    names(fl) <- nm
  } else {
#    warning('The fl-argument is obsolete')
  }

  if(!is.list(fl)) stop('need at least one factor')
  fl <- lapply(fl,as.factor)
  if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')

#  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf,'terms')

  y <- model.response(mf,'numeric')

# try a sparse model matrix to save memory when removing intercept
# though, demeanlist must be full.  Ah, no, not much to save because
# it won't be sparse after centering
# we should rather let demeanlist remove the intercept, this
# will save memory by not copying.  But we need to remove it below in x %*% beta
# (or should we extend beta with a zero at the right place, it's only
#  a vector, eh, is it, do we not allow matrix lhs? No.)

  x <- model.matrix(mt,mf)
  rm(mf)

  icpt <- 0
  icpt <- which(attr(x,'assign') == 0)
  if(length(icpt) == 0) icpt <- 0
  ncov <- ncol(x) - (icpt > 0)
  if(ncov == 0) {
    # No covariates
    fr <- demeanlist(y,fl)
    z <- list(r.residuals=y,fe=fl,p=0,cfactor=compfactor(fl),residuals=fr,call=match.call())
    class(z) <- 'felm'
    return(z)
  }
  # here we need to demean things

  dm <- demeanlist(list(y=y,x=x),fl,icpt)
  yz <- dm[[1]]
  xz <- dm[[2]]
  rm(dm)
  gc()

  badconv <- attr(xz,'badconv') + attr(yz,'badconv')
  dim(xz) <- c(nrow(x),ncov)
  attributes(yz) <- attributes(y)


# here we just do an lm.fit, however lm.fit is quite slow since
# it doesn't use blas (in particular it can't use e.g. threaded blas in acml)
# so we have rolled our own.

# we really don't return an 'lm' object or other similar stuff, so
# we should consider using more elementary operations which map to blas-3
# eg. solve(crossprod(xz),t(xz) %*% yz)
# Or, even invert by solve(crossprod(xz)) since we need
# the diagonal for standard errors.  We could use the cholesky inversion
# chol2inv(chol(crossprod(xz)))

  cp <- crossprod(xz)
  ch <- cholx(cp)

#  ch <- chol(cp)
#  beta <- drop(inv %*% (t(xz) %*% yz))
  # remove multicollinearities
  badvars <- attr(ch,'badvars')
  b <- crossprod(xz,yz)

  if(is.null(badvars)) {
    beta <- as.vector(backsolve(ch,backsolve(ch,b,transpose=TRUE)))
    inv <- chol2inv(ch)
  } else {
    beta <- rep(NaN,nrow(cp))
    beta[-badvars] <- backsolve(ch,backsolve(ch,b[-badvars],transpose=TRUE))
    inv <- matrix(NaN,nrow(cp),ncol(cp))
    inv[-badvars,-badvars] <- chol2inv(ch)
  }
  rm(b)
  if(icpt > 0) names(beta) <- colnames(x)[-icpt] else names(beta) <- colnames(x)
#  cat(date(),'projected system finished\n')
  z <- list(coefficients=beta,badconv=badconv)
  N <- nrow(xz)
  p <- ncol(xz) - length(badvars)

# how well would we fit with all the dummies?
# the residuals of the centered model equals the residuals
# of the full model, thus we may compute the fitted values
# resulting from the full model.
  zfit <- xz %*% ifelse(is.na(beta),0,beta)

  rm(xz)
  zresid <- yz - zfit
  rm(yz)
  z$fitted.values <- y - zresid
  z$residuals <- zresid
  # insert a zero at the intercept position
  if(length(fl) > 0) {
    if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta
    pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
    z$r.residuals <- y - pred
  } else {
    z$r.residuals <- zresid
  }
#  z$xb <- pred
  rm(x)
  rm(y)
  gc()

  z$cfactor <- compfactor(fl)

  numrefs <- nlevels(z$cfactor) + max(length(fl)-2,0)
  numdum <- sum(unlist(lapply(fl,nlevels))) - numrefs
  z$numrefs <- numrefs
#  if(length(fl) <= 2) {
#    numdum <- sum(unlist(lapply(fl,nlevels))) - nlevels(z$cfactor)
#  } else {
#    numdum <- sum(unlist(lapply(fl,nlevels))) - length(fl) + 1
#  }
  z$df <- N - p - numdum
  vcvfactor <- sum(z$residuals**2)/z$df
  z$vcv <- inv * vcvfactor
  z$se <- sqrt(diag(z$vcv))
  z$sefactor <- sqrt(vcvfactor)
  z$tval <- z$coefficients/z$se
  z$pval <- 2*pt(abs(z$tval),z$df,lower.tail=FALSE)
  z$terms <- mt
  z$fe <- fl
  z$N <- N
  z$p <- p + numdum
  z$xp <- p
  z$call <- match.call()
  class(z) <- 'felm'

  return(z)
}
