# $Id: oldfelm.R 1942 2016-04-07 21:25:18Z sgaure $
# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0




#   Some things in this file is done in a weird way.
#   In some cases there are efficiency reasons for this, e.g. because
#   the "standard" way of doing things may result in a copy which is costly
#   when the problem is *large*.
#   In other cases it may simply be due to the author's unfamiliarity with how
#   things should be done in R
 
# parse our formula
oldparseformula <- function(formula,data) {

  trm <- terms(formula,specials=c('G'))
  feidx <- attr(trm,'specials')$G+1
  va <- attr(trm,'variables')
  festr <- paste(sapply(feidx,function(i) deparse(va[[i]])),collapse='+')

  if(festr != '') { 
    .Deprecated(msg="The G() syntax is deprecated, please use multipart formulas instead")
    # remove the G-terms from formula
    formula <- update(formula,paste('. ~ . -(',festr,') - 1'))
    
    # then make a list of them, and find their names
    felist <- parse(text=paste('list(',gsub('+',',',festr,fixed=TRUE),')',sep=''))
    nm <- eval(felist,list(G=function(arg) deparse(substitute(arg))))
    
    # replace G with factor, eval with this, and the parent frame, or with data
    # allow interaction factors with '*' (dropped, never documented, use ':')
    Gfunc <- function(f) if(is.null(attr(f,'xnam'))) factor(f) else f
    Ginfunc <- function(x,f) {
      if(is.factor(x)) {
        structure(interaction(factor(f),factor(x),drop=TRUE),xnam=deparse(substitute(x)),fnam=deparse(substitute(f)))
      } else {
        structure(factor(f),x=x,xnam=deparse(substitute(x)), fnam=deparse(substitute(f)))
      }
    }
    
    if(is.environment(data)) {
      fl <- eval(felist,list(G=Gfunc, ':'=Ginfunc),data)
    } else {
      fl <- local({eval(felist,data)},list(G=Gfunc, ':'=Ginfunc))
    }
    names(fl) <- nm
    gpart <- eval(parse(text=paste('~',paste(nm,collapse='+'))))

    if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')
  } else {
    fl <- NULL
    gpart <- ~0
  }
  return(list(formula=formula, fl=fl, gpart=gpart,ivpart=~0,cpart=~0))
}

# parse 
# use 2-part Formulas without G() syntax, like
# y ~ x1 + x2 | f1+f2
# or 3-part or more Formulas with iv-specification like
# y ~ x1 + x2 | f1+f2 | (q|w ~ x3+x4) | c1+c2
# returns a list containing
# formula=y~x1+x2
# fl = list(f1,f2)
# ivpart = list(q ~x3+x4, w ~x3+x4)
# cluster=list(c1,c2)
nopart <- function(x) length(all.vars(x))==0

parseformula <- function(form, data, noexpand=FALSE) {
  f <- as.Formula(form)
  len <- length(f)[[2]]
  if(len == 1) return(oldparseformula(form,data))
  opart <- formula(f,lhs=NULL,rhs=1)
  if(len == 1) return(list(formula=opart,gpart=~0,ivpart=~0,cpart=~0))

  # the factor part
  gpart <- formula(f,lhs=0, rhs=2)
  if(!nopart(gpart)) {
    tm <- terms(gpart, keep.order=TRUE)
    ft <- attr(tm,'factors')
    var <- eval(attr(tm,'variables'),data)
    varnames <- rownames(ft)
    names(var) <- varnames
    fl <- apply(ft, 2, function(v) {
      nonz <- sum(v > 0)
      vnam <- varnames[which(v > 0)]
      if(nonz > 2) stop('Interaction only supported for two variables')
      if(nonz == 1) {
#        if(!is.factor(var[[vnam]])) warning('non-factor ',vnam, ' coerced to factor')
        res <- list(factor(var[[vnam]]))
        names(res) <- vnam
      } else {
        xnam <- vnam[[1]]
        fnam <- vnam[[2]]
        x <- var[[xnam]]
        f <- var[[fnam]]
        if(!is.factor(f) && !is.factor(x)) {
          stop('interaction between ', xnam, ' and ', fnam, ', none of which are factors')
        }
        if(!is.factor(f) && is.factor(x)) {
          tmp <- x
          x <- f
          f <- tmp
          tmp <- xnam
          xnam <- fnam
          fnam <- tmp
        }
        if(is.factor(x)) {
          res <- list(structure(interaction(factor(f),factor(x),drop=TRUE),xnam=xnam,fnam=fnam))
        } else {
          res <- list(structure(factor(f),x=x,xnam=xnam, fnam=fnam))
        }
        names(res) <- paste(xnam,fnam,sep=':')
      }
      res
    })
    nm <- names(fl)
    fl <- unlist(fl, recursive=FALSE)
    names(fl) <- nm
  } else {
    fl <- NULL
  }
  
  if(len == 2) return(list(formula=opart,fl=fl,gpart=gpart,ivpart=~0,cpart=~0))

  # Then the iv-part
  ivparts <- formula(f,lhs=0,rhs=3, drop=TRUE)
  if(!nopart(ivparts) && length(ivparts[[2]])>1 && ivparts[[2]][[1]]=='(') {
      # Now, make a list of the iv-formulas where we split the lhs in each
      # to obtain q ~ x3+x4, w ~x3+x4
      ivspec <- as.Formula(ivparts[[2]][[2]])   # it's now q|w ~ x3+x4
      lhs <- formula(ivspec,rhs=0)
      ivpart <- lapply(seq_along(all.vars(lhs)),function(i) formula(ivspec,lhs=i))
  } else {
      ivpart <- NULL
  }

  if(len == 3 && !is.null(ivpart)) return(list(formula=opart,fl=fl,iv=ivpart,gpart=gpart,ivpart=ivparts,cpart=~0))
  
  # The cluster part, this could be the third part if there are no parentheses
  if(len == 3 && is.null(ivpart)) {
    cpart <- ivparts
    ivparts <- NULL
  } else {
    cpart <- formula(f,lhs=0,rhs=4,drop=TRUE)
  }
  if(!nopart(cpart)) {
      # handle the same way as the factors, but without the covariate interaction
      tm <- terms(cpart, keep.order=TRUE)
      nm <- parts <- attr(tm,'term.labels')
      clist <- lapply(paste('factor(',parts,')',sep=''),function(e) parse(text=e))
      cluster <- lapply(clist,eval,data)
      names(cluster) <- nm
      
  } else {
      cluster <- NULL
  }
  list(formula=opart,fl=fl,iv=ivpart,cluster=cluster,gpart=gpart,ivpart=ivparts,cpart=cpart)
}


# ivresid is optional, used in 2. stage of 2sls to pass
# the difference between the original endogenous variable and the prediction
# for the purpose of computing sum of square residuals
doprojols <- function(psys, ivresid=NULL, exactDOF=FALSE, keepX=FALSE, nostats=FALSE) {
  if(is.numeric(exactDOF)) {
    df <- exactDOF
    totvar <- length(psys$y) - df
  } else {
    # numrefs is also used later
    numrefs <- nrefs(psys$fl, compfactor(psys$fl), exactDOF) 
    totvar <- totalpvar(psys$fl)-numrefs
    df <- length(psys$y)-totvar
  }
  if(is.null(psys$yxz$x)) {
    # No covariates
        z <- list(N=psys$N, r.residuals=psys$y,fe=psys$fl,p=totvar,Pp=0,cfactor=compfactor(psys$fl),
                  na.action=psys$na.action, contrasts=psys$contrasts,
                  fitted.values=psys$y - psys$yxz$y,
                  coefficients=matrix(double(0),psys$N,0),
                  df=df,
                  nostats=FALSE,
                  model.assign=psys$assign,
                  model.labels=psys$model.labels,
                  residuals=psys$yxz$y,clustervar=psys$clustervar, call=match.call())
        z$df.residual <- z$df
        class(z) <- 'felm'
        return(z)
    }
  yz <- psys$yxz$y
  xz <- psys$yxz$x
  y <- psys$y
  x <- psys$x
  fl <- psys$fl
  icpt <- psys$icpt
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
  b <- crossprod(xz,yz)
  ch <- cholx(cp)
    #  ch <- chol(cp)
    #  beta <- drop(inv %*% (t(xz) %*% yz))
    # remove multicollinearities
  badvars <- attr(ch,'badvars')

  if(is.null(badvars)) {
    beta <- backsolve(ch,backsolve(ch,b,transpose=TRUE))
#    beta <- as.vector(beta)
#    beta <- as.vector(backsolve(ch,backsolve(ch,b,transpose=TRUE)))
    if(!nostats) inv <- chol2inv(ch)
  } else {
    beta <- matrix(NaN, nrow(cp), ncol(b))
#    beta <- rep(NaN,nrow(cp))
    beta[-badvars,] <- backsolve(ch,backsolve(ch,b[-badvars,], transpose=TRUE))
    if(!nostats) {
      inv <- matrix(NA,nrow(cp),ncol(cp))
      inv[-badvars,-badvars] <- chol2inv(ch)
    }
  }
  rm(ch, b, cp)

    
  if(length(fl) > 0 && icpt > 0) 
    rownames(beta) <- colnames(x)[-icpt] else rownames(beta) <- colnames(x)
  colnames(beta) <- colnames(y)
  if(ncol(beta) == 1) names(beta) <- rownames(beta)
  
  z <- list(coefficients=beta,badconv=psys$badconv,Pp=ncol(xz))


  z$N <- nrow(xz)
  z$p <- ncol(xz) - length(badvars)
  if(!nostats) {
    z$inv <- inv
    inv <- nazero(inv)
  }
# how well would we fit with all the dummies?
# the residuals of the centered model equals the residuals
# of the full model, thus we may compute the fitted values
# resulting from the full model.

# for the 2. step in the 2sls, we should replace 
# the instrumented variable with the real ones (the difference is in ivresid)
# when predicting, but only for the purpose of computing
# residuals.

  nabeta <- nazero(beta)

  zfit <- xz %*% nabeta
  zresid <- yz - zfit
  z$beta <- beta
  z$response <- y
  z$fitted.values <- y - zresid
  z$residuals <- zresid
  z$contrasts <- psys$contrasts
  z$model.assign <- psys$model.assign
  z$model.labels <- psys$model.labels
  if(length(fl) > 0) {
  # insert a zero at the intercept position (x may have an intercept, whereas xz has not)
#    if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta
    if(icpt > 0) {
      pre <- seq_len(icpt-1)
      post <- setdiff(seq_len(nrow(beta)), pre)
      ibeta <- rbind(beta[pre,,drop=FALSE], 0, beta[post,,drop=FALSE])
    } else {
      ibeta <- beta
    }

    pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
    z$r.residuals <- y - pred
  } else {
    z$r.residuals <- zresid
  }
  rm(x)


  z$lhs <- colnames(beta)

  # the residuals should be the residuals from the original endogenous variables, not the predicted ones
  # the difference are the ivresid, which we must multiply by beta and subtract.
  # the residuals from the 2nd stage are in iv.residuals
  # hmm, what about the r.residuals?  We modify them as well. They are used in kaczmarz().
  if(!is.null(ivresid)) {
    if(!is.matrix(ivresid)) {
      nm <- names(ivresid)
      ivresid <- matrix(unlist(ivresid),z$N)
      colnames(ivresid) <- nm
    }
    z$ivresid <- ivresid %*% nabeta[colnames(ivresid),,drop=FALSE]
    z$iv.residuals <- z$residuals
    z$residuals <- z$residuals - z$ivresid    
    z$r.iv.residuals <- z$r.residuals
    z$r.residuals <- z$r.residuals - z$ivresid
  }
  
  z$terms <- psys$terms
  z$cfactor <- compfactor(fl)
  totlev <- totalpvar(fl)
  if(is.numeric(exactDOF)) {
    z$df <- exactDOF
    numdum <- z$N - z$p - z$df
    z$numrefs <- totlev - numdum
  } else {
    numdum <- totlev - numrefs
    z$numrefs <- numrefs
    z$df <- z$N - z$p - numdum
  }
  z$df.residual <- z$df
  z$rank <- z$N - z$df
  
  z$exactDOF <- exactDOF

  z$fe <- fl
# should we subtract 1 for an intercept?
# a similar adjustment is done in summary.felm when computing rdf
  z$p <- z$p + numdum - 1  
  z$xp <- z$p
  z$na.action <- psys$na.action
  class(z) <- 'felm'
  cluster <- psys$clustervar
  z$clustervar <- cluster

  if(nostats) {
    z$nostats <- TRUE
    return(z)
  }

  z$nostats <- FALSE
# then we go about creating the covariance matrices and tests
# if there is a single lhs, they are just stored as matrices etc
# in z.  If there are multiple lhs, these quantities are inserted
# in a list z$STATS indexed by z$lhs
# indexed by the name of the lhs

  vcvnames <- list(rownames(beta), rownames(beta))
  Ncoef <- nrow(beta)

  singlelhs <- length(z$lhs) == 1
  if(!singlelhs) z$STATS <- list()

  for(lhs in z$lhs) {
    res <- z$residuals[,lhs]

#  if(!is.null(ivresid)) res <- res - z$ivresid

    vcvfactor <- sum(res**2)/z$df

# when multiple lhs, vcvfactor is a vector
# we need a list of vcvs in this case

    if(singlelhs) {
      z$vcv <- z$inv * vcvfactor
      setdimnames(z$vcv, vcvnames)
    } else {
      z$STATS[[lhs]] <- list()
      z$STATS[[lhs]]$vcv <- z$inv * vcvfactor
      setdimnames(z$STATS[[lhs]]$vcv, vcvnames)
    }

#  dimnames(z$vcv) <- list(names(beta),names(beta))

  # We should make the robust covariance matrix too.
  # it's inv * sum (X_i' u_i u_i' X_i) * inv
  # where u_i are the (full) residuals (Wooldridge, 10.5.4 (10.59))
  # i.e. inv * sum(u_i^2 X_i' X_i) * inv
  # for large datasets the sum is probably best computed by a series of scaled
  # rank k updates, i.e. the dsyrk blas routine, we make an R-version of it.
  # need to check this computation, the SE's are slightly numerically different from Stata's.
  # it seems stata does not do the small-sample adjustment
    dfadj <- z$N/z$df

  # Now, here's an optimzation for very large xz. If we use the wcrossprod and ccrossprod
  # functions, we can't get rid of xz, we end up with a copy of it which blows away memory.
  # we need to scale xz with the residuals in xz, but we don't want to expand res to a full matrix,
  # and even get a copy in the result.
  # Thus we modify it in place with a .Call. The scaled variant is also used in the cluster computation.
#  z$robustvcv <- dfadj * inv %*% wcrossprod(xz,res) %*% inv

    rscale <- ifelse(res==0,1e-40,res)
    .Call(C_scalecols, xz, rscale)
    if(singlelhs) {
      z$robustvcv <- dfadj * inv %*% crossprod(xz) %*% inv
      setdimnames(z$robustvcv, vcvnames)
    } else {
      z$STATS[[lhs]]$robustvcv <- dfadj * inv %*% crossprod(xz) %*% inv
      setdimnames(z$STATS[[lhs]]$robustvcv, vcvnames)
    }


  # then the clustered covariance matrix
    if(!is.null(cluster)) {
      method <- attr(cluster,'method')
      if(is.null(method)) method <- 'cgm'
      dfadj <- (z$N-1)/z$df
      d <- length(cluster)
      if(method == 'cgm') {
#        meat <- matrix(0,nrow(z$vcv),ncol(z$vcv))
        meat <- matrix(0,Ncoef,Ncoef)
        for(i in 1:(2^d-1)) {
          # Find out which ones to interact
          iac <- as.logical(intToBits(i))[1:d]
          # odd number is positive, even is negative
          sgn <- 2*(sum(iac) %% 2) - 1
          # interact the factors
          ia <- factor(do.call(paste,c(cluster[iac],sep='\004')))
          adj <- sgn*dfadj*nlevels(ia)/(nlevels(ia)-1)
          .Call(C_dsyrk,1,meat,adj,rowsum(xz,ia))
        }
        if(singlelhs) {
          z$clustervcv <- inv %*% meat %*% inv
          setdimnames(z$clustervcv, vcvnames)
        } else {
          z$STATS[[lhs]]$clustervcv <- inv %*% meat %*% inv
          setdimnames(z$STATS[[lhs]]$clustervcv, vcvnames)
        }
        rm(meat)
        ## } else if(method == 'gaure') {
        ##   .Call(C_scalecols, xz, 1/rscale)
        ##   meat <- matrix(0,nrow(z$vcv),ncol(z$vcv))
        ##   dm.res <- demeanlist(res,cluster)
        ##   skel <- lapply(cluster, function(f) rep(0,nlevels(f)))
        ##   means <- relist(kaczmarz(cluster,res-dm.res), skel)
        ##   scale <- ifelse(dm.res==0,1e-40, dm.res)
        ##   .Call(C_scalecols, xz, scale)
        ##   .Call(C_dsyrk, 1, meat, dfadj, xz)
        ##   .Call(C_scalecols, xz, 1/scale)
        ##   for(i in seq_along(cluster)) {
        ##     rs <- rowsum(xz, cluster[[i]])
        ##     adj <- nlevels(cluster[[i]])/(nlevels(cluster[[i]])-1)
        ##     .Call(C_scalecols, rs, means[[i]])
        ##     .Call(C_dsyrk, 1, meat, dfadj*adj, rs)
        ##   }
        ##   rm(xz,rs)
        ##   z$clustervcv <- inv %*% meat %*% inv
        ##   rm(meat)
      } else {
        stop('unknown multi way cluster algorithm:',method)
      }
      
      
      if(singlelhs) {
        z$cse <- sqrt(diag(z$clustervcv))
        z$ctval <- coef(z)/z$cse
        z$cpval <- 2*pt(abs(z$ctval),z$df,lower.tail=FALSE)
      } else {
        z$STATS[[lhs]]$cse <- sqrt(diag(z$STATS[[lhs]]$clustervcv))
        z$STATS[[lhs]]$ctval <- z$coefficients[,lhs]/z$STATS[[lhs]]$cse
        z$STATS[[lhs]]$cpval <- 2*pt(abs(z$STATS[[lhs]]$ctval),z$df,lower.tail=FALSE)
      }
    }
    if(singlelhs) {
      z$se <- sqrt(diag(z$vcv))
      z$tval <- z$coefficients/z$se
      z$pval <- 2*pt(abs(z$tval),z$df,lower.tail=FALSE)

      z$rse <- sqrt(diag(z$robustvcv))
      z$rtval <- coef(z)/z$rse
      z$rpval <- 2*pt(abs(z$rtval),z$df,lower.tail=FALSE)
    } else {
      z$STATS[[lhs]]$se <- sqrt(diag(z$STATS[[lhs]]$vcv))
      z$STATS[[lhs]]$tval <- z$coefficients[,lhs]/z$STATS[[lhs]]$se
      z$STATS[[lhs]]$pval <- 2*pt(abs(z$STATS[[lhs]]$tval),z$df,lower.tail=FALSE)

      z$STATS[[lhs]]$rse <- sqrt(diag(z$STATS[[lhs]]$robustvcv))
      z$STATS[[lhs]]$rtval <- z$coefficients[,lhs]/z$STATS[[lhs]]$rse
      z$STATS[[lhs]]$rpval <- 2*pt(abs(z$STATS[[lhs]]$rtval),z$df,lower.tail=FALSE)
    }
    # reset this for next lhs
    if(!singlelhs) .Call(C_scalecols, xz, 1/rscale)
  }
  
  z
}


project <- function(mf,fl,data,contrasts,clustervar=NULL,pf=parent.frame()) {

  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  subspec <- mf[['subset']]

  # we should handle multiple lhs
  # but how?  model.frame() doesn't handle it, but we need
  # model.frame for subsetting and na.action, with the left hand side
  # included.  We create an artifical single lhs by summing the left hand
  # sides, just to get hold of the rhs.  Then we extract the left hand side

  Form <- as.Formula(mf[['formula']])
  mf[['formula']] <- Form
  mf <- eval(mf, pf)
  mt <- attr(mf,'terms')
  naact <- attr(mf,'na.action')
  if(!is.null(naact))
    naclass <- class(naact)
  fullN <- nrow(mf) + length(naact)
  # then obtain the response matrix through Formula::model.part
  response <- as.matrix(model.part(Form, mf, lhs=NULL,rhs=0))


  cmethod <- attr(clustervar,'method')
  if(!is.null(clustervar)) {
      if(is.character(clustervar)) clustervar <- as.list(clustervar)
      if(!is.list(clustervar)) clustervar <- list(clustervar)
      clustervar <- lapply(clustervar, function(cv) {
          if(!is.character(cv)) factor(cv) else factor(data[,cv])
      })
  }
  # we need to change clustervar and factor list to reflect
  # subsetting and na.action. na.action is ok, it's set as an attribute in mf
  # but subset must be done manually. It's done before na handling
  if(!is.null(subspec)) {
    subs <- eval(subspec,pf)
    if(!is.null(clustervar)) clustervar <- lapply(clustervar,function(cv) cv[subs])
    fl <- lapply(fl, function(fac) {
      f <- factor(fac[subs])
      x <- attr(f,'x')
      if(is.null(x)) return(f)
      structure(f,x=x[subs])
    })
                 
  }
  if(!is.null(naact)) {
    if(!is.null(clustervar)) clustervar <- lapply(clustervar, function(cv) cv[-naact])
    fl <- lapply(fl,function(fac) {
      f <- factor(fac[-naact])
      x <- attr(f,'x')
      if(is.null(x)) return(f)
      structure(f,x=x[-naact])
    })
  }

  attr(clustervar,'method') <- cmethod

#  ret <- list(fl=fl, na.action=naact,terms=mt,clustervar=clustervar, y=model.response(mf,'numeric'))
  ret <- list(fl=fl, na.action=naact, fullN=fullN, terms=mt,clustervar=clustervar, y=response)  
  rm(mt,clustervar,naact)

  lapply(ret$clustervar, function(f)
         if(length(f) != nrow(ret$y)) stop('cluster factors are not the same length as data ',
                      length(f),'!=',nrow(ret$y)))

  # in case of cluster factor specified with the clustervar argument:

# try a sparse model matrix to save memory when removing intercept
# though, demeanlist must be full.  Ah, no, not much to save because
# it won't be sparse after centering
# we should rather let demeanlist remove the intercept, this
# will save memory by not copying.  But we need to remove it below in x %*% beta
# (or should we extend beta with a zero at the right place, it's only
#  a vector, eh, is it, do we not allow matrix lhs? No.)

# we make some effort to avoid copying the data matrix below
# this includes assigning to lists in steps, with gc() here and there.
# It's done for R 3.0.2. The copy semantics could be changed in later versions.

#  ret$x <- model.matrix(ret$terms,mf,contrasts)
  ret$x <- model.matrix(Form, mf, contrasts)
  rm(mf)
  ret$contrasts <- attr(ret$x,'contrasts')
  ret$model.assign <- attr(ret$x, 'assign')
  ret$model.labels <- attr(terms(Form[-2]), 'term.labels')   # ditch lhs when finding terms
  icpt <- attr(ret$x,'assign') == 0
  if(!any(icpt)) icpt <- 0 else icpt <- which(icpt)
  ret$icpt <- icpt

  ncov <- ncol(ret$x) - (icpt > 0)
  if(ncov == 0) {
    ret$x <- NULL
    ret$yxz <- list(y=demeanlist(ret$y,fl))
    ret$Pp <- 0
    ret$N <- length(ret$y)
    ret$yx <- list(y=ret$y)
    return(ret)
  }

  # here we need to demean things
  # we take some care so that unexpected copies don't occur
  # hmm, the list() copies the stuff. How can we avoid a copy
  # and still enable parallelization over y and x in demeanlist? A vararg demeanlist?
  # I.e. an .External version?  Yes.
#  yx <- list(y=ret$y, x=ret$x)
#  gc()
#  ret$yxz <- demeanlist(yx,fl,icpt)
#  rm(fl,yx); gc()

#  ret$yxz <- edemeanlist(y=ret$y,x=ret$x,fl=fl,icpt=c(0,icpt))
  ret$yxz <- demeanlist(list(y=ret$y,x=ret$x),fl=fl,icpt=c(0,icpt))
  ret$badconv <- attr(ret$yxz$x,'badconv') + attr(ret$yxz$y,'badconv')
  # use our homebrewn setdimnames instead of colnames. colnames copies.
  if(length(fl) > 0) {
    if(icpt == 0)
      setdimnames(ret$yxz$x, list(NULL,colnames(ret$x)))
    else 
      setdimnames(ret$yxz$x, list(NULL,colnames(ret$x)[-icpt]))
  }

  ret
}

#' @export
..oldfelm <- function(formula, data, exactDOF=FALSE, subset, na.action, contrasts=NULL,...) {

  knownargs <- c('iv', 'clustervar', 'cmethod', 'keepX', 'nostats')
  keepX <- FALSE
  cmethod <- 'cgm'
  iv <- NULL
  clustervar <- NULL
  nostats <- FALSE

  deprec <- c('iv', 'clustervar')

#  sc <- names(sys.call())[-1]
#  named <- knownargs[pmatch(sc,knownargs)]
#  for(arg in c('iv', 'clustervar')) {
#    if(!is.null(eval(as.name(arg))) && !(arg %in% named)) {
#        warning("Please specify the '",arg,"' argument by name, or use a multi part formula. Its position in the argument list will change in a later version")
#      }
#  }

  mf <- match.call(expand.dots = TRUE)

  # Currently there shouldn't be any ... arguments
  # check that the list is empty

#  if(length(mf[['...']]) > 0) stop('unknown argument ',mf['...'])
  
  # When moved to the ... list, we use this:
  # we do it right away, iv and clustervar can't possibly end up in ... yet, not with normal users


  args <- list(...)
  ka <- knownargs[pmatch(names(args),knownargs, duplicates.ok=FALSE)]
  names(args)[!is.na(ka)] <- ka[!is.na(ka)]
  dpr <- deprec[match(ka, deprec)]
  if(any(!is.na(dpr))) {
    bad <- dpr[which(!is.na(dpr))]
    warning('Argument(s) ',paste(bad,collapse=','), ' are deprecated and will be removed, use multipart formula instead')
  }
  env <- environment()
  lapply(intersect(knownargs,ka), function(arg) assign(arg,args[[arg]], pos=env))
  if(!(cmethod %in% c('cgm','gaure'))) stop('Unknown cmethod: ',cmethod)

  # also implement a check for unknown arguments
  unk <- setdiff(names(args), knownargs)
  if(length(unk) > 0) stop('unknown arguments ',paste(unk, collapse=' '))

  if(missing(data)) mf$data <- data <- environment(formula)
  pf <- parent.frame()
  pform <- parseformula(formula,data)
  
  if(!is.null(iv) && !is.null(pform[['iv']])) stop("Specify EITHER iv argument(deprecated) OR multipart terms, not both")
  if(!is.null(pform[['cluster']]) && !is.null(clustervar)) stop("Specify EITHER clustervar(deprecated) OR multipart terms, not both")
  if(!is.null(pform[['cluster']])) clustervar <- structure(pform[['cluster']], method=cmethod)
  if(is.null(iv) && is.null(pform[['iv']])) {
    # no iv, just do the thing
    fl <- pform[['fl']]
    formula <- pform[['formula']]
    mf[['formula']] <- formula
    psys <- project(mf,fl,data,contrasts,clustervar,pf)
    z <- doprojols(psys,exactDOF=exactDOF, nostats=nostats[1])
    if(keepX) z$X <- if(psys$icpt > 0) psys$x[,-psys$icpt] else psys$x
    rm(psys)

    z$parent.frame <- pf
    z$call <- match.call()
    return(z)
  }

  # IV.  Clean up formulas, set up for 1st stages
  if(!is.null(iv)) {
    # warning("argument iv is deprecated, use multipart formula instead")
    if(!is.list(iv)) iv <- list(iv)
    form <- pform[['formula']]
    # Old syntax, the IV-variables are also in the main equation, remove them
    for(ivv in iv) {
      ivnam <- ivv[[2]]
      # create the new formula by removing the IV lhs.
      form <- update(form, substitute(. ~ . - Z,list(Z=ivnam)))
    }
    pform[['formula']] <- form
    mf[['iv']] <- NULL
    ivpart <- NULL
  } else {
    iv <- pform[['iv']]
    ivpart <- as.Formula(pform[['ivpart']])
  }


  if(is.environment(data)) {
    ivenv <- new.env(parent=data)
  } else {
    ivenv <- new.env(parent=pf)
  }


  if(!is.null(ivpart)) {
    # ivpart is something like ~(P|Q ~ x+x2)
    # strip ~ and ()
    ivpart <- as.Formula(ivpart[[2]][[2]])
    lhs <- formula(ivpart, lhs=NULL, rhs=0)
    rhs <- as.Formula(formula(ivpart, lhs=0, rhs=1))
    if(length(ivpart)[[2]] > 1) {
      stop("Instruments can't be projected out: ", ivpart)
      rhsg <- as.Formula(formula(ivpart, lhs=0, rhs=2))
    } else {
      rhsg <- list(NULL,NULL)
    }
    Form <- as.Formula(formula)
    if(length(Form)[2] == 4) {
      cluform <- formula(Form,lhs=0,rhs=4)[[2]]
      step1form <- as.Formula(substitute(L ~ B + ivO | G|0|C,
                                         list(L=lhs[[2]], B=pform[['formula']][[3]],
                                              ivO=rhs[[2]],
                                              G=pform[['gpart']][[2]],
                                              C=cluform)))
      
    } else {
      step1form <- as.Formula(substitute(L ~ B + ivO | G,
                                         list(L=lhs[[2]], B=pform[['formula']][[3]],
                                              ivO=rhs[[2]],
                                              G=pform[['gpart']][[2]])))
    }
    environment(step1form) <- environment(formula)


    ivform <- parseformula(step1form,data)
    fl <- ivform[['fl']]
    mf[['formula']] <- ivform[['formula']]
    environment(mf[['formula']]) <- environment(formula)
    psys <- project(mf,fl,data,contrasts,clustervar,pf)

    if(length(nostats) == 2)
        nost <- nostats[2]
    else
        nost <- nostats[1]
    z1 <- doprojols(psys, exactDOF=exactDOF, nostats=nost)
    # put the fitted values in ivenv
    for(n in colnames(z1$fitted.values)) {
      nn <- paste(n,'(fit)',sep='')
      assign(nn,z1$fitted.values[,n],envir=ivenv)
    }
    z1$endogvars <- paste('`',colnames(z1$fitted.values),'(fit)`',sep='')

    # we should not use all.vars(rhs), to find the instruments, but
    # pick up things from z1 somehow, in case there are expanded
    # factors as instrumental variables
    # in z1 there are model.assign and model.labels.
    # we locate the instrument names (rhs) in model.labels, and extract the coefficient names from model.assign
    cname <- rownames(z1$coefficients)
    ivnam <- all.vars(rhs)
    asgn <- z1$model.assign
    if(length(asgn) != length(cname)) asgn <- asgn[asgn!=0]
    # now assume there's a factor f among the instruments, with levels f2-f4 (1 is a reference)
    # we look up 'f' in lab, it has position 5, then everything with assign == 5 belong to this factor
    ivars <- cname[asgn %in% which(z1$model.labels %in% ivnam)]
    z1$instruments <- ivars
    # Store any dummy instruments in the ivenv, for possible later use in condfstat
    # There's a psys$x which is the model matrix
    # we may as well store all the instruments there
    for(ivar in ivars) {
      if(!(ivar %in% ivnam))
          assign(ivar,psys$x[,ivar],envir=ivenv)
    }
    if(!nost) {
      z1$iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1,ivars, lhs=lh))
      names(z1$iv1fstat) <- z1$lhs
      z1$rob.iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1,ivars,type='robust', lhs=lh))
      names(z1$rob.iv1fstat) <- z1$lhs
    }

    z1$call <- match.call()
    environment(step1form) <- environment(formula)
    z1$call[['formula']] <- step1form

    naact <- psys$na.action

    FIT <- as.formula(paste('~',paste(z1$endogvars,sep='', collapse='+'),sep=''))
    step2form <- as.Formula(substitute(y ~ B + FIT | G,
                                    list(y=pform[['formula']][[2]],
                                         B=pform[['formula']][[3]],
                                         FIT=FIT[[2]],
                                         G=pform[['gpart']][[2]])))

    # do the first step
    environment(step2form) <- environment(formula)

    form2 <- parseformula(step2form, data)
    fl <- form2[['fl']]
    formula <- form2[['formula']]
    environment(formula) <- ivenv
    mf$formula <- formula
    if(is.environment(mf$data)) mf$data <- ivenv
    # remove the naact from 1st step
    # any missings in the exogenous variables will be missing in both the 1st and 2nd stage
    # If there are missing instruments or endogenous variables they will be missing in 1st stage, but not in the 2nd
    # so we must make them missing in the 2nd stage as well. We implement that as a subset.
    if(!is.null(naact)) {
      # naact is a vector of observations to remove
      # numbered after a subset has been taken.
      # if there's no subset in mf, we create one
      subset <- mf[['subset']]
      if(is.null(subset)) {
        mf[['subset']] <- -naact
      } else {
        # subset is an indexing vector for the rows in the data frame
        # naact specifies rows to remove after subsetting.
        # we simulate this process. Let's make an integer vector for the
        # rows of the data. The total length can be found
        mf[['subset']] <- seq_len(psys$fullN)[subset][-naact]
      }
    }
    psys <- project(mf=mf,fl=fl,data=data,contrasts=contrasts,clustervar=clustervar,pf=pf)
    ivres <- z1$residuals
    colnames(ivres) <- as.character(sapply(all.vars(FIT),as.name))
    z <- doprojols(psys, ivresid=ivres, exactDOF=exactDOF, nostats=nostats)
    z$stage1 <- z1
    z$st2call <- mf
# backwards compatibility
    z$step1 <- lapply(z1$lhs, function(lh) {
#      warning('Use stage1 instead of step1')
      foo <- z1
      foo$lhs <- lh
      foo$beta <- z1$beta[,lh,drop=FALSE]
      foo$coefficients <- foo$beta
      foo$response <- z1$response[,lh,drop=FALSE]
      foo$fitted.values <- z1$fitted.values[,lh,drop=FALSE]
      foo$residuals <- z1$residuals[,lh,drop=FALSE]
      foo$r.residuals <- z1$r.residuals[,lh,drop=FALSE]
      if(!is.null(z1$iv1fstat)) {
        foo$iv1fstat <- z1$iv1fstat[lh]
        foo$rob.iv1fstat <- z1$rob.iv1fstat[lh]
      }
      if(!is.null(z1$STATS)) 
          foo[names(z1$STATS[[lh]])] <- z1$STATS[[lh]]
      foo[['STATS']] <- NULL
      foo})
    z$endovars <- z1$endogvars
    z$parent.frame <- pf
    rm(psys)
    z$call <- match.call()
    return(z)
  }

  # parse the IV-formulas, they may contain factor-parts
  iv <- lapply(iv,parseformula,data)

  # Now, insert the rhs of the IV in the formula
  # find the ordinary and the factor part

  # It may contain a factor part

  # this is the template for the step1 formula, just insert a left hand side
  step1form <- formula(as.Formula(substitute(~B +ivO | G + ivG,
                                    list(B=pform[['formula']][[3]],G=pform[['gpart']][[2]]))))

  # this is the template for the second stage, it is updated with the IV variables
  step2form <- formula(as.Formula(substitute(y~B | G,
                                    list(y=pform[['formula']][[2]],B=pform[['formula']][[3]],
                                         G=pform[['gpart']][[2]]))))

  nullbase <- formula(as.Formula(substitute(~B | G,
                                    list(B=pform[['formula']][[3]],G=pform[['gpart']][[2]]))))

  # we must do the 1. step for each instrumented variable
  # collect the instrumented variables and remove them from origform
  # then do the sequence of 1. steps
  # we put the instrumented variable on the lhs, and add in the formula for it on the rhs
  # A problem with this approach is that na.action is run independently between the first stages and
  # the second stage. This may result in different number of observations. This isn't any good.
  # I haven't figured out a good solution for this, except for creating the full model.matrix for
  # all of the steps.  This will blow our memory on large datasets.

  ivarg <- list()
  vars <- NULL
  step1 <- list()
  endolist <- c()
  for(ivv in iv) {
    # Now, make the full instrumental formula, i.e. with the rhs expanded with the
    # instruments, and the lhs equal to the instrumented variable

     ivlhs <- ivv[['formula']][[2]]
     rhsivo <- formula(as.Formula(ivv[['formula']]),lhs=0)[[2]]
    if(nopart(ivv[['gpart']]))
        rhsivg <- 0
    else
        rhsivg <- formula(ivv[['gpart']],lhs=0,rhs=2)[[2]]

    fformula <- substitute(Z ~ R, list(Z=ivlhs, R=step1form[[2]]))

    fformula <- do.call(substitute, list(fformula,list(ivO=rhsivo,ivG=rhsivg)))

     ivform <- parseformula(fformula,data)
     fl <- ivform[['fl']]
     mf[['formula']] <- ivform[['formula']]

     # note that if there are no G() terms among the instrument variables,
     # all the other covariates should only be centered once, not in every first stage and the
     # second stage separately. We should rewrite and optimize for this.
     psys <- project(mf,fl,data,contrasts,clustervar,pf)
     z <- doprojols(psys, exactDOF=exactDOF)
     mf[['formula']] <- fformula
     z$call <- mf
     rm(psys)


     # now, we need an ftest between the first step with and without the instruments(null-model)
     # We need the residuals with and without the
     # instruments. We have them with the instruments, but must do another estimation
     # for the null-model
     mfnull <- mf
     nullform <- substitute(Z ~ R,list(Z=ivlhs,R=nullbase[[2]]))     
     pformnull <- parseformula(nullform, data)
     mfnull[['formula']] <- pformnull[['formula']]
     znull <- doprojols(project(mfnull, pformnull[['fl']], data, contrasts, clustervar, pf),
                        exactDOF=exactDOF)
     z$iv1fstat <- ftest(z,znull)
     z$rob.iv1fstat <- ftest(z,znull,vcov=z$robustvcv)
     if(!is.null(clustervar))
       z$clu.iv1fstat <- ftest(z,znull,vcov=z$clustervcv)
     step1 <- c(step1,list(z))
     # then we lift the fitted variable and create a new name
     ivz <- z
     evar <- deparse(ivlhs)
     new.var <- paste(evar,'(fit)',sep='')
     # store them in an environment
     assign(new.var,ivz$fitted.values,envir=ivenv)
#     data[[new.var]] <- ivz$fitted
     # save these, with the backtick for later use
     vars <- c(vars,paste('`',new.var,'`',sep=''))
     # keep the residuals, they are needed to reconstruct the residuals for the
     # original variables in the 2. stage
     ivarg[[paste('`',new.var,'`',sep='')]] <- ivz$residuals
     # and add it to the equation
     step2form <- update(as.Formula(step2form),as.formula(substitute(. ~ . + FIT | .,
                                         list(FIT=as.name(new.var)))))
  }
  names(step1) <- names(iv)
  # now we have a formula in step2form with all the iv-variables
  # it's just to project it

  pform <- parseformula(step2form,data)
  fl <- pform[['fl']]
  formula <- pform[['formula']]
  environment(formula) <- ivenv
  mf$formula <- formula
  if(is.environment(mf$data)) mf$data <- ivenv
  psys <- project(mf=mf,fl=fl,data=data,contrasts=contrasts,clustervar=clustervar,pf=pf)

  z <- doprojols(psys,ivresid=ivarg,exactDOF=exactDOF)
  z$step1 <- step1
  z$endovars <- vars
  rm(psys)
  z$call <- match.call()
  return(z)
}



