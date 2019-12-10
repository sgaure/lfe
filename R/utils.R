# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0

# this call changes or adds the dimnames of obj without duplicating it
# it should only be used on objects with a single reference
# Our use of it is safe, and we don't export it.
setdimnames <- function(obj, nm) {
  .Call(C_setdimnames,obj,nm)
}

scalecols <- function(obj, vec) {
  .Call(C_scalecols, obj, vec)
}

Crowsum <- function(x,f,mean=FALSE) .Call(C_rowsum,x,f,mean)

orthonormalize <- function(V) {
  structure(V %*% solve(chol(crossprod(V))), ortho=TRUE)
}

# This function sets an attribute 'inplace' on the argument
# It is checked in the C-function demeanlist, and then removed
# If set, the centering will be inplace. It can mess things
# up seriously if used in the wrong manner, so it's not public
# except as in an argument to demeanlist (through eval-trickery).
unnamed <- function(x) {
  .Call(C_inplace,x)
}

wmean <- function(x,w) sum(w*x)/sum(w)
wcov <- function(x,y,w) {
  sw <- w/sum(w)
  mx <- sum(sw*x)
  my <- sum(sw*y)
  cx <- x - mx
  cy <- y - my
  sum(sw*cx*cy)/(1-sum(sw^2))
}
wvar <- function(x,w) wcov(x,x,w)

pinvx <- function(X) {
#    return(pinv(nazero(X)))
    ch <- cholx(nazero(X))
    badvars <- attr(ch,'badvars')
    inv1 <- chol2inv(ch)
    if(is.null(badvars)) return(inv1)
    inv <- matrix(0,nrow(X),ncol(X))
    inv[-badvars,-badvars] <- inv1
    structure(inv,badvars=attr(ch,'badvars'))
}

# Do a Cholesky to detect multi-collinearities
cholx <- function(mat, eps=1e-6) {
  if(is.null(dim(mat))) dim(mat) <- c(1,1)
  N <- dim(mat)[1]
  if(N == 1) {
      return(structure(sqrt(mat),badvars=if(mat<=0) 1 else NULL))
  }

  # first, try a pivoted one
  tol <- N*getOption('lfe.eps')
  chp <- chol(mat,pivot=TRUE,tol=tol)
  rank <- attr(chp,'rank')
  if(rank == N) return(chol(mat))
  pivot <- attr(chp,'pivot')
  oo <- order(pivot)
  badvars <- pivot[((rank+1):N)]
  ok <- (1:N)[-badvars]
  ch <- chol(mat[ok,ok])
  return(structure(ch,badvars=badvars))
}


cholsolve <- function(A,b) {
  ch <- chol(A)
  backsolve(ch,backsolve(ch,b,transpose=TRUE))
}

pinv <- function (X, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  badvars <- integer(0)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  Positive <- ifelse(is.na(Positive),FALSE,Positive)
  badvars <- which(!Positive)
  if (all(Positive))
    res <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    res <- array(0, dim(X)[2L:1L])
  else
    res <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                            t(Xsvd$u[, Positive, drop = FALSE]))
  structure(res,badvars=badvars)
}

nazero <- function(x) ifelse(is.na(x),0,x)

ftest <- function(z, zr=NULL, vcov=z$vcv) {
  rdf <- z$df #z$N - z$p - 1
  if(is.null(zr)) {
# we should do F-test vs. model with intercept.
# but the intercept in z is implicit.
    F <- (t(coef(z)) %*% pinvx(vcov) %*% coef(z))/z$p
    return( c(F=F, p=pf(F, z$p, rdf, lower.tail=FALSE), df1=z$p, df2=rdf))
  }
#  df1 <- length(coef(z)) - length(coef(zr))
  df1 <- z$p - zr$p
  c1 <- coef(z)
  c2 <- rep(0,length(c1))
  names(c2) <- names(c1)
  c2[names(coef(zr))] <- coef(zr)
  F <- (t(c1-c2) %*% pinvx(vcov) %*% (c1-c2))/df1
  return(c(F=F, p=pf(F,df1, rdf, lower.tail=FALSE), df1=df1,df2=rdf))
}


mkgraph <- function(f1,f2) {
  if(!requireNamespace('igraph', quietly=TRUE)) stop('Package igraph not found')
#  graph.edgelist(cbind(paste('f1',f1),paste('f2',f2)), directed=FALSE)
#  graph.edgelist(cbind(500000000+as.integer(f1),f2), directed=FALSE)
  igraph::graph.adjacency(tcrossprod(makeDmatrix(list(f1,f2)))>0,
                  'undirected', diag=FALSE)
}

diamgraph <- function(flist,approx=TRUE) {
  if(!requireNamespace('igraph', quietly=TRUE)) stop('Package igraph not found')
  gr <- mkgraph(flist[[1]],flist[[2]])
# find largest cluster
  cl <- igraph::clusters(gr)$membership
  lcl <- which(cl == which.max(table(cl)))
  if(approx) 
    max(igraph::shortest.paths(gr,v=sample(lcl,10),to=sample(lcl,10)))
  else
    igraph::diameter(igraph::induced.subgraph(gr,lcl))
}







#' Find diameters of mobility graphs
#' 
#' 'diammatrix' computes the diameters of certain graphs related to convergence
#' speed of \code{felm}.
#' 
#' Each pair of factors (f1,f2) from \code{flist} defines a bipartite graph in
#' which the vertices are the levels of the factors, and two vertices are
#' adjacent if they are observed simultaneously. The connected components of
#' this graph are important for identification of the coefficients for the
#' factor levels, i.e. for \code{getfe}. But experience and some trials have
#' led the author to speculate that the diameter of the graph (or its largest
#' component) is also important for the convergence rate.  Specifically, the
#' author suspects that under some assumptions, time to convergence goes like
#' the square of the diameter.  At least in the case of two factors.  This
#' function computes the diameter for each pair of factors.  If the graph is
#' disconnected, the largest connected component is used. If \code{accel=TRUE}
#' (the default), the diameter is approximated from below by drawing two sets
#' of 10 random vertices and finding the maximum length of the shortest paths
#' between them.
#' 
#' @param flist a list of factors defining the dummies.
#' @param approx logical. Approximate diameters are computed.
#' @return A matrix of dimension K x K where K is \code{length(flist)}.
#' @note This function is not important to the operation of the package, it is
#' included for easy experimentation with the convergence rate.  It requires
#' that the suggested package \pkg{igraph} is attached.
#' @keywords internal
diammatrix <- function(flist, approx=TRUE) {

  flen <- length(flist)
  if(flen < 2) return(0)
  val <- matrix(0,flen,flen)
  colnames(val) <- names(flist)
  rownames(val) <- names(flist)
  for(if1 in 1:(flen-1))
    for(if2 in (if1+1):flen) {
      val[if1,if2] <- diamgraph(flist[c(if1,if2)], approx)
    }

# fill in lower triangle:
  val[row(val) > col(val)] <- t(val)[row(val) > col(val)]
  val
}

# if(!exists('fixef')) fixef <- function(object,...) UseMethod('fixef')

# compute rank deficiency of D-matrix
rankDefic <- function(fl,method='cholesky',mctol=1e-3) {
  eps <- sqrt(.Machine$double.eps)
  if(length(fl) == 1) return(1)
  if(length(fl) == 2) return(nlevels(compfactor(fl)))
  if(method == 'cholesky') {
    D <- makeDmatrix(fl)
    Ch <- as(Cholesky(crossprod(D), super=TRUE, perm=TRUE, Imult=eps), 'sparseMatrix')
    sum(diag(Ch) < eps^(1/4))
  } else if(method == 'mc') {
    totlev <- sum(sapply(fl,nlevels))
    N <- length(fl[[1]])
    len <- length(fl) + 1
    # now, we will use the rank deficiency d to compute
    # degrees of freedom corrections. I.e. (tr+totlev)+len should
    # be within a certain relative tolerance, which translates
    # to an absolute tolerance of tr
    tolfun <- function(tr) -mctol*abs((tr+totlev)+len)
    init <- 0
    N - (mctrace(fl,tol=tolfun,init=init) + totlev)+len
    
  } else {
    D <- makeDmatrix(fl)
    as.integer(ncol(D) - rankMatrix(crossprod(D), method='qr.R'))
  }
}


# in makeDmatrix/makePmatrix, the weights argument should be the
# square root of the weights. This is what is stored in internal structures.
# 






#' Make sparse matrix of dummies from factor list
#' 
#' Given a list of factors, return the matrix of dummies as a sparse matrix.
#' 
#' The function returns the model matrix for a list of factors. This matrix is
#' not used internally by the package, but it's used in some of the
#' documentation for illustrative purposes.
#' 
#' @param fl list of factors.
#' @param weights numeric vector. Multiplied into the rows.
#' @return Returns a sparse matrix.
#' @examples
#' 
#'   fl <- lapply(1:3, function(i) factor(sample(3,10,replace=TRUE)))
#'   fl
#'   makeDmatrix(fl, weights=seq(0.1,1,0.1))
#' 
#' @export makeDmatrix
makeDmatrix <- function(fl, weights=NULL) {
  # make the D matrix
  # for pure factors f, it's just the t(as(f,'sparseMatrix'))
  # if there's a covariate vector x, it's t(as(f,'sparseMatrix'))*x
  # for a covariate matrix x, it's the cbinds of the columns with the factor-matrix
  ans <- do.call(mycbind,lapply(fl, function(f) {
    x <- attr(f,'x',exact=TRUE)
    fm <- t(as(f,'sparseMatrix'))
    if(is.null(x)) return(fm)
    if(!is.matrix(x)) return(fm*x)
    do.call(mycbind,apply(x,2,'*',fm))
  }))
  nm <- names(fl)
  if(is.null(nm)) nm <- paste('f',seq_along(fl),sep='')
  levnm <- unlist(sapply(seq_along(fl), function(i) xlevels(nm[i],fl[[i]])))
  if(!is.null(weights)) 
      ans <- Diagonal(length(weights),weights) %*% ans
  colnames(ans) <- levnm
  ans
}

makePmatrix <- function(fl, weights=NULL) {
  D <- makeDmatrix(fl,weights)
  DtD <- crossprod(D)
  DtDi <- pinvx(as.matrix(DtD))
  badvars <- attr(DtDi,'badvars')
  if(is.null(badvars)) return(D %*% DtDi %*% t(D))
  D <- D[,-badvars,drop=FALSE]
  return(D %*% pinvx(as.matrix(crossprod(D))) %*% t(D))
}
# total number of variables projected out
totalpvar <- function(fl) {
  if(length(fl) == 0) return(0)
  sum(sapply(fl, function(f) {
    x <- attr(f,'x',exact=TRUE)
    if(is.null(x) || !is.matrix(x)) return(nlevels(f))
    return(ncol(x)*nlevels(f))
  }))
}

nrefs <- function(fl, cf, exactDOF=FALSE) {
  if(length(fl) == 0) return(0)
  if(missing(cf)) cf <- compfactor(fl)
  numpure <- sum(sapply(fl,function(f) is.null(attr(f,'x',exact=TRUE))))
  if(numpure == 1) return(0)
  if(numpure == 2) return(nlevels(cf))
  if(identical(exactDOF,'rM')) {
    return(rankDefic(fl, method='qr'))
  } else if(identical(exactDOF,'mc')) {
    return(rankDefic(fl, method='mc'))
  } else  if(exactDOF) {
    return(rankDefic(fl, method='cholesky'))
  }
  return(nlevels(cf) + numpure-2)
}

wildcard <- function(formula,  s=ls(environment(formula)), re=FALSE,
                     nomatch.failure=TRUE) {
  env <- environment(formula)
  F <- as.Formula(formula)
  lenlhs <- length(F)[1]
  lenrhs <- length(F)[2]
  if(lenlhs == 1 && lenrhs == 1) return(swildcard(formula,s,re, nomatch.failure=nomatch.failure));
  # do the parts separately
  lhs <- lapply(seq_len(lenlhs), function(lh) {
    swildcard(formula(F, lhs=lh,rhs=0), s, re, nomatch.failure=nomatch.failure)[[2]]
  })
  rhs <- lapply(seq_len(lenrhs), function(rh) {
    swildcard(formula(F, lhs=0,rhs=rh), s, re, nomatch.failure=nomatch.failure)[[2]]
  })
  if(length(lhs) > 0) {
    LHS <- lhs[[1]]
    for(s in lhs[-1]) {
      LHS <- bquote(.(LHS) | .(s))
    }
  } else LHS <- NULL

  RHS <- rhs[[1]]
  for(s in rhs[-1]) {
    RHS <- bquote(.(RHS) | .(s))
  }
  if(is.null(LHS)) return(as.Formula(substitute( ~ R, list(R=RHS)), env=env))
   else
  as.Formula(substitute(L ~ R, list(L=LHS, R=RHS)), env=env)
}

swildcard <- function(formula, s=ls(environment(formula)), re=FALSE,
                      nomatch.failure=TRUE) {
  av <- all.vars(formula)
  if(!re) {
    rchars <- c('*','?')
    wild <- unique(unlist(sapply(rchars,grep, av,fixed=TRUE,value=TRUE)))
    if(length(wild) > 0)
        rewild <- utils::glob2rx(wild,trim.tail=FALSE)
    else
        rewild <- NULL
    # quote some regexp-specials
    rewild <- lapply(seq_along(rewild), function(i) structure(rewild[i],orig=wild[i]))
  } else {
    rchars <- c('.','*','|','\\','?','[',']','(',')','{','}')
    wild <- unique(unlist(sapply(rchars, grep, av, fixed=TRUE, value=TRUE)))
    rewild <- lapply(wild, function(w) structure(w,orig=w))
  }
  wsub <- lapply(rewild, function(w) {
    mtch <- grep(paste('^',w,'$',sep=''), s, value=TRUE, perl=TRUE)
    if(length(mtch) == 0) {
      if(nomatch.failure) stop("Couldn't match wildcard variable `",
                               attr(w,'orig'),'`', call.=FALSE)
      mtch <- attr(w,'orig')
    }
    as.list(parse(text=paste('`',mtch,'`',collapse='+',sep='')))[[1]]
  })
  names(wsub) <- wild
  formula(as.Formula(do.call(substitute, list(formula,wsub)),
                     env=environment(formula)),
          update=T, collapse=TRUE, drop=FALSE)
}


#prettyprint a list of integers
ilpretty <- function(il) {
  paste(tapply(il,cumsum(c(1,diff(il)) != 1),
         function(x) if(length(x) == 1) {
           as.character(x)
         } else {
           paste(x[[1]],x[[length(x)]],sep=':')
         }), collapse=' ')
}

makefitnames <- function(s) gsub('`(Intercept)(fit)`','(Intercept)',paste('`',s,'(fit)`',sep=''), fixed=TRUE)
delete.icpt <- function(x) {
  asgn <- attr(x,'assign')
  icpt <- asgn == 0
  if(!length(icpt)) return(x)
  asgn <- asgn[!icpt]
  ctr <- attr(x,'contrasts')
  x <- x[,!icpt,drop=FALSE]
  attr(x,'assign') <- asgn
  attr(x,'contrasts') <- ctr
  x
}

# do an ls on env, and the parent and all the way up to top
rls <- function(env, top=.GlobalEnv) {
  ret <- character()
  e <- env
  repeat {
    if(identical(e,emptyenv())) break;
    ret <- c(ret,ls(e))
    if(identical(e,top)) break;
    e <- parent.env(e)
  }
  unique(ret)
}

limlk <- function(mm) {
  #  find the kappa for computing k-class liml
  # it's the smallest eigenvalue of the matrix
  # [y endog]' M_i [y endog] ([y endog]' M [y endog])^{-1}
  # where M is projection out of instruments
  # and M_i is projection out of both instruments and exogenous
  # variables
  ye <- cbind(mm$y, mm$ivy)
  H1 <- crossprod(ye, newols(list(y=ye,x=mm$x),nostats=TRUE)$residuals)
  H <- crossprod(ye, newols(list(y=ye,x=cbind(mm$ivx,mm$x)), nostats=TRUE)$residuals)
  mat <- solve(H,H1)
  min(eigen(mat, only.values=TRUE)$values)
}

resample <- function(cluster, return.factor=FALSE, na.action=NULL, fill=FALSE) {
  if(is.list(cluster)) {
    if(is.factor(cluster[[1]])) {
      if(length(cluster) > 1) warning('List of factors not supported in resample, using first only (',names(cluster)[1],')')
      cluster <- cluster[[1]]
    }
  }
  # resample entire levels of a factor
  if(length(na.action) > 0) cluster <- factor(cluster[-na.action])
  iclu <- as.integer(cluster)
  cl <- sort(sample(nlevels(cluster), replace=TRUE))
  if(fill) {
    tb <- table(cluster)
    while(sum(tb[cl]) < length(cluster)) {
      cl <- c(cl,sample(nlevels(cluster),1))
    }
    cl <- sort(cl)
  }
  # find a faster way to do this:
  # s <- sort(unlist(sapply(cl, function(ll) which(clu==ll))))
  s <- NULL
  while(length(cl) > 0) {
    s <- c(s,which(iclu %in% cl))
    cl <- cl[c(1L,diff(cl)) == 0]
  }
  if(return.factor) return(cluster[s])
  sort(s)
}

getcomp <- function(est, alpha=NULL) {
  if(nlevels(est$cfactor) == 1) return(list(est=est,alpha=alpha))
  ok <- which(est$cfactor == 1)
  res <- est
  res$cfactor <- factor(res$cfactor[ok])
  if(!is.null(res$X))
      res$X <- res$X[ok,]
  res$residuals <- res$residuals[ok,]
  res$response <- res$response[ok,]
  res$fitted.values <- res$fitted.values[ok,]
  res$r.residuals <- res$r.residuals[ok,]
  res$fe <- lapply(est$fe, function(f) factor(f[ok]))
  if(alpha != NULL) {
    # of the two first factors, remove the components beyond the first
    # if there are more than two factors, the remaining are assumed to be ok
    # so find those with 'fe' among the two first, and comp > 1
    
    bad <- which((alpha[,'fe'] %in% names(est$fe)[1:2]) & (alpha[,'comp'] > 1))
    if(length(bad) > 0) alpha <- alpha[-bad,]
  }
  return(est=res,alpha=alpha)
}

#' Chain subset conditions 
#'
#' @param ... Logical conditions to be chained. 
#' @param out.vars character. Variables not in data.frame, only needed if you use variables which
#' are not in the frame.  If \code{out.vars} is not specified, it is assumed to match all variables
#' starting with a dot ('.').
#' @return Expression that can be \code{eval}'ed to yield a logical subset mask.
#' @details
#' A set of logical conditions are chained, not and'ed. That is, each argument to
#' \code{chainsubset} is used as a filter to create a smaller dataset. Each subsequent
#' argument filters further.
#' For independent conditions this will be the same as and'ing them. I.e.
#' \code{chainsubset(x < 0 , y < 0)} will yield  the same subset as \code{(x < 0) & (y < 0)}.
#' However, for aggregate filters like \code{chainsubset(x < mean(y), x > mean(y))}
#' we first find all the observations with \code{x < mean(y)}, then among these we
#' find the ones with \code{x > mean(y)}.  The last \code{mean(y)} is now conditional on
#' \code{x < mean(y)}.
#' 
#' @examples
#' set.seed(48)
#' N <- 10000
#' dat <- data.frame(y=rnorm(N), x=rnorm(N))
#' # It's not the same as and'ing the conditions:
#' felm(y ~ x,data=dat,subset=chainsubset(x < mean(y), y < 2*mean(x)))
#' felm(y ~ x,data=dat,subset=chainsubset(y < 2*mean(x), x < mean(y)))
#' felm(y ~ x,data=dat,subset=(x < mean(y)) & (y < 2*mean(x)))
#' lm(y ~ x, data=dat, subset=chainsubset(x < mean(y), x > mean(y)))
#' @note
#' Some trickery is done to make this work directly in the subset argument of functions like
#' \code{felm()} and \code{lm()}. It might possibly fail with an error message in some situations.
#' If this happens, it should be done in two steps: \code{ss <- eval(chainsubset(...),data); 
#' lm(...,data=data, subset=ss)}. In particular, the arguments are taken literally, 
#' constructions like \code{function(...) {chainsubset(...)}} or \code{a <- quote(x < y); chainsubset(a)} do
#' not work, but \code{do.call(chainsubset,list(a))} does.
#' @export
chainsubset <- function(..., out.vars) {
  if(sys.parent() != 0) {
    cl <- sys.call(sys.parent())
    mc <- match.call(expand.dots=TRUE)
    if(cl[[1]] == quote(eval)) {
      # we are called from eval, probably inside lm. let's replace it with a evalq and a direct call, then
      # evaluation of the result
      ecaller <- parent.frame(3)
      cl[[1]] <- quote(evalq)
      cl[[2]] <- match.call()
      cl[[2]] <- eval(cl,ecaller)
      return(eval(cl,ecaller))
    }
  }

  args <- as.list(match.call(expand.dots=TRUE))[-1]
  args[['out.vars']] <- NULL
  if(length(args) < 1) return(NULL)
  if(length(args) == 1) return(structure(args[[1]],filter=args))
  group <- list(as.name('{'))
  miss <- missing(out.vars)
  R <- as.name(sprintf('R%x',RR <- sample(.Machine$integer.max,1)))
  FF <- as.name(sprintf('F%x',RR))
  structure(
      bquote(local({
        .(R) <- logical(length(.(FF) <- .(args[[1]])))
        .(R)[.(Reduce(function(f1,f2) {
          nm <- all.vars(f2)
          nm <- lapply(if(miss)
                         grep('^\\.',nm,value=TRUE,invert=TRUE)
                       else
                         nm[!(nm %in% out.vars)],as.name)
          recod <- as.call(c(group, lapply(nm, function(n) bquote(.(n) <- .(n)[.(FF), drop=TRUE]))))
          bquote({.(FF) <- .(f1); local({.(recod); .(FF)[.(f2)]})})
        }, args[-1],init=bquote(base::which(.(FF)))))] <- TRUE
        .(R)
      })),
      filter=args)
}
