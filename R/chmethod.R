# $Id: chmethod.R 1693 2015-04-07 09:36:29Z sgaure $
findfe <- function(dd,Rhs,se=FALSE) {
  # find references
  refnames <- attr(dd,'refnames')
  nm <- c(colnames(dd),refnames)
  refcnt <- attr(dd,'refcnt')
  # add in the reference in dg

  dg <- c(diag(dd),refcnt)
  ok <- 1:ncol(dd)

  if(se) sev <- double(length(nm))
  alphacoef <- double(length(nm))

  # the super-nodal algorithm
  # is default and far better, but it consumes more memory

  trysolve <- try(solve(dd,Rhs))
  if(inherits(trysolve,'try-error')) {
    if(grepl('problem too large',geterrmessage())) {
      message(paste('Never mind, trying *slower* non-supernodal algorithm, nnz=',nnzero(dd)))
      message(paste(date(),'This may be an opportunity for a nice cup of tea. Or two.'))
      gc()
      ch <- Cholesky(dd,super=FALSE,perm=TRUE)
      trysolve <- solve(ch,Rhs)
      rm(ch)
      gc()
    } else {
      stop(geterrmessage())
    }
  } 


  alphacoef[ok] <- as.vector(trysolve)
  if(se) {
    # is there a faster way to find the diagonal of the inverse?
    sev[ok] <- sqrt(diag(solve(dd)))*attr(se,'sefactor')
    alpha <- data.frame(effect=alphacoef,se=sev,obs=dg)
  } else {
    alpha <- data.frame(effect=alphacoef,obs=dg)
  }
  rownames(alpha) <- nm
  alpha
}


#makedummies <- function(factors) {
#  nm <- c()
#  dummies <- Matrix(0,0,length(factors[[1]]))
#  for(i in 1:length(factors)) {
#    f <- factors[[i]]
#    dummies <- rBind(dummies,as(f,'sparseMatrix'))
#    nm <- c(nm,paste(names(factors)[[i]],levels(f),sep='.'))
#  }
#  rownames(dummies) <- nm
#  dummies
#}

makedd.full <- function(factors) {
#  dm <- makedummies(factors)
  dm <- t(makeDmatrix(factors))
  nm <- rownames(dm)
  dd <- tcrossprod(dm)
  rownames(dd) <- colnames(dd) <- nm
  attr(dd,'dummies') <- dm
  attr(dd,'nm') <- nm
  dd
}

makeddlist <- function(factors) {
  if(length(factors) > 2) {
    if(is.null(attr(factors,'references'))) {
      # find references by fiddling with Cholesky
      message('*** More than two groups, finding refs by Cholesky pivots, interpret at own risk')
      # first the full matrix, find small pivots
      dd <- makedd.full(factors)
      orignm <- attr(dd,'nm')

      # add small amount to diagonal
      eps <- sqrt(.Machine$double.eps)
      Ch <- try(Cholesky(dd,super=TRUE,perm=TRUE,Imult=eps))
      if(inherits(Ch,'try-error') && grepl('problem too large',geterrmessage())) {
        Ch <- Cholesky(dd,super=FALSE,perm=TRUE,Imult=eps)
      }
      # strangely enough, coercing to sparseMatrix doesn't take care of
      # the permutation, we apply it manually.  Let's hope it's never fixed.
      rm(dd); gc()
      pivot <- Ch@perm
      ch <- as(Ch,'sparseMatrix')
      rm(Ch); gc()
      dg <- diag(ch)[order(pivot)]**2
      rm(ch); gc()
      refs <- (dg < eps**(1/3))
      refnames <- orignm[refs]
      message(paste('***',length(refnames),'references found'))
    } else {
      refnames <- attr(factors,'references')
      orignm <- unlist(lapply(names(factors),
                              function(n) paste(n,levels(factors[[n]]),sep='.')))
    }
    # there may be references in more than one factor
    # remove all of them
    # create factor list with named levels
    nf <- lapply(names(factors),function(n) {
      f <- factors[[n]]
      levels(f) <- paste(n,levels(f),sep='.')
      f
    })
    # remove reference levels, and remove the prefix
    # find the levels
    lev <- lapply(nf,function(f) which(levels(f) %in% refnames))
    nnf <- mapply(function(f,l) factor(f,exclude=levels(f)[l]),factors,lev,SIMPLIFY=FALSE)
    dd <- makedd.full(nnf)
    attr(dd,'keep') <- 1:length(nnf[[1]])
    attr(dd,'refnames') <- refnames
#    attr(dd,'refcnt') <- rep(1,length(refnames))
# find the number of occurences
    cntlst <- unlist(lapply(refnames,function(n) lapply(nf,function(f) sum(f == n))))
    attr(dd,'refcnt') <- cntlst[cntlst > 0]
    attr(dd,'comp') <- 1
    res <- list(dd)
    attr(res,'nm') <- orignm
  } else {
    # 2 or fewer factors, find references by component
    cf <- compfactor(factors)
    nml <- lapply(factors,function(f) levels(f))
    nm <- unlist(lapply(names(nml),function(n) paste(n,nml[[n]],sep='.')))
    res <- list()
    li <- 1
# this loop suffers from too much copying and stuff
# when there are many components (e.g. like 10000)
    remfact <- factors
    fullidx <- 1:length(factors[[1]])
    for(l in levels(cf)) {
      # find those in this level
      keep <- which(cf == l)
#      cat(date(),'comp fact',li,'size',length(keep),'\n')
      fcomp <- lapply(remfact,function(f) factor(f[keep]))
      remfact <- lapply(remfact,function(f) factor(f[-keep]))
      cf <- factor(cf[-keep])
      # then the reference level
      maxrefs <- lapply(fcomp,function(f) {tf <- table(f); m <- which.max(tf); tf[m]})
      # in which factor
      rfac <- which.max(unlist(maxrefs))
      # which level
      reflevel <- names(maxrefs[[rfac]])
      # drop that level from the factor
      fcomp[[rfac]] <- factor(fcomp[[rfac]],exclude=reflevel)
      refname <- paste(names(remfact)[[rfac]],reflevel,sep='.')
      # remove those without levels
      len <- unlist(lapply(fcomp,nlevels))
      fcomp <- fcomp[len > 0]
      dd <- makedd.full(fcomp)
      # the keep attribute should be relative to
      # the full factor, not to remfact
      attr(dd,'keep') <- fullidx[keep]
      fullidx <- fullidx[-keep]
      attr(dd,'refnames') <- refname
      attr(dd,'refcnt') <- max(unlist(maxrefs))
      attr(dd,'comp') <- li
      res[[li]] <- dd
      li <- li+1
#      res <- c(res, list(dd))
    }
    attr(res,'nm') <- nm
  }
  res
}
