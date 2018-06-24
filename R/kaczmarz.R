#' Solve a linear system defined by factors
#' 
#' Uses the Kaczmarz method to solve a system of the type Dx = R, where D is
#' the matrix of dummies created from a list of factors.
#' 
#' 
#' @param fl A list of arbitrary factors of the same length
#' @param R numeric.  A vector, matrix or list of such of the same length as
#' the factors
#' @param eps a tolerance for the method
#' @param init numeric. A vector to use as initial value for the Kaczmarz
#' iterations. The algorithm converges to the solution closest to this
#' @param threads integer. The number of threads to use when \code{R} is more
#' than one vector
#' @return A vector \code{x} of length equal to the sum of the number of levels
#' of the factors in \code{fl}, which solves the system \eqn{Dx=R}. If the
#' system is inconsistent, the algorithm may not converge, it will give a
#' warning and return something which may or may not be close to a solution. By
#' setting \code{eps=0}, maximum accuracy (with convergence warning) will be
#' achieved.
#' @note This function is used by \code{\link{getfe}}, it's quite specialized,
#' but it might be useful for other purposes too.
#' 
#' In case of convergence problems, setting \code{options(lfe.usecg=TRUE)} will
#' cause the kaczmarz() function to dispatch to the more general conjugate
#' gradient method of \code{\link{cgsolve}}.  This may or may not be faster.
#' @seealso \code{\link{cgsolve}}
#' @examples
#' 
#' ## create factors
#'   f1 <- factor(sample(24000,100000,replace=TRUE))
#'   f2 <- factor(sample(20000,length(f1),replace=TRUE))
#'   f3 <- factor(sample(10000,length(f1),replace=TRUE))
#'   f4 <- factor(sample(8000,length(f1),replace=TRUE))
#' ## the matrix of dummies
#'   D <- makeDmatrix(list(f1,f2,f3,f4))
#'   dim(D)
#' ## an x
#'   truex <- runif(ncol(D))
#' ## and the right hand side
#'   R <- as.vector(D %*% truex)
#' ## solve it
#'   sol <- kaczmarz(list(f1,f2,f3,f4),R)
#' ## verify that the solution solves the system Dx = R
#'   sqrt(sum((D %*% sol - R)^2))
#' ## but the solution is not equal to the true x, because the system is
#' ## underdetermined
#'   sqrt(sum((sol - truex)^2))
#' ## moreover, the solution from kaczmarz has smaller norm
#'   sqrt(sum(sol^2)) < sqrt(sum(truex^2))
#' 
#' @export kaczmarz
kaczmarz <- function(fl,R,eps=getOption('lfe.eps'),init=NULL,
                     threads=getOption('lfe.threads')) {

  if(getOption('lfe.usecg')) {
    mat <- makeDmatrix(fl)
    if(is.list(R)) {
      mm <- crossprod(mat)
      return(lapply(R, function(ll) drop(cgsolve(mm, crossprod(mat, ll), eps=max(eps,1e-6), init=init))))
#      skel <- lapply(R,function(a) {if(is.matrix(a)) matrix(0,ncol(mat),ncol(a)) else rep(0,ncol(mat))})
#      return(utils::relist(cgsolve(crossprod(mat), crossprod(mat,Reduce(cbind,R)), eps=eps, init=init), skel))
    }
    return(drop(cgsolve(crossprod(mat), crossprod(mat,R), eps=eps, init=init)))
  }
  if(is.null(threads)) threads <- 1
  islist <- is.list(R)
  if(!islist) R <- list(R)
  v <- .Call(C_kaczmarz,fl,R,eps,as.vector(init),as.integer(threads))
  if(!islist) {
    v <- drop(v[[1]])
  }
  v
}


getfe.kaczmarz <- function(obj,se=FALSE,eps=getOption('lfe.eps'),ef='ref',bN=100,
                           robust=FALSE, cluster=NULL, lhs=NULL) {
  if(is.character(ef)) {
    ef <- efactory(obj,opt=ef)
  }
  if(!isTRUE(attr(ef,'verified')) && !is.estimable(ef, obj$fe)) {
    warning('Supplied function seems non-estimable')
  }
  multlhs <- length(obj$lhs) > 1
  if(is.null(lhs)) {
    R <- obj$r.residuals-obj$residuals
  } else {
    if(!all(lhs %in% obj$lhs)) stop('lhs must be subset of ', paste(obj$lhs, collapse=' '))
    R <- obj$r.residuals[,lhs, drop=FALSE] - obj$residuals[,lhs, drop=FALSE]
  }
  v <- kaczmarz(obj$fe,R,eps)

  if(is.matrix(v) && ncol(v) > 1) {
    v <- apply(v,2,ef,addnames=TRUE)
    vtmp <- ef(v[,1],addnames=TRUE)
    extra <- attr(vtmp, 'extra')
    nm <- names(vtmp)
  } else {
    v <- ef(v,TRUE)
    extra <- attr(v,'extra')
    nm <- names(v)
  }

  res <- data.frame(effect=v)
  if(multlhs) colnames(res) <- paste('effect',colnames(R),sep='.')
  if(!is.null(extra)) res <- cbind(res,extra)
  rownames(res) <- nm
  attr(res,'ef') <- ef

  if(se) {
    if(multlhs) {
      for(lh in colnames(R)) {
        res <- btrap(res,obj,bN,eps=eps, robust=robust, cluster=cluster, lhs=lh)        
      }
    } else {
      res <- btrap(res,obj,bN,eps=eps, robust=robust, cluster=cluster)
    }
  }
  res
}


# A common estimable function on the fe-coefficients
# return an estimable function, the matrix which
# pins a reference in each component, the one with the
# most observations
# if there are more than two factors, assume they don't
# ruin identification beyond one extra reference for each such factor

# Note: when I get some time, I'll implement a Weeks-Williams estimable
# function.  I.e. with a reference in each factor in each component
# from compfactor(...,WW=TRUE). This may be what many people want.
# No, can't do that. The same level may occur in separate WW-components.
# WW is about partitioning of the dataset, not of the levels.

# return level names in appropriate order
# the G(x:f) with x a matrix makes it slightly complicated
xlevels <- function(n,f,sep='.') {
  x <- attr(f,'x',exact=TRUE)
  plev <- paste(n,levels(f),sep=sep)
  if(is.null(x) || !is.matrix(x)) return(plev)
  nam <- attr(f,'xnam')
  if(is.null(nam)) nam <- 'x'
  if(!is.matrix(x)) return(paste(nam,plev,sep=sep))
  matnam <- colnames(x)
  if(is.null(matnam)) matnam <- paste(nam,1:ncol(x),sep='') else matnam <- paste(nam,matnam,sep='')
  plev <- sub('.*:','',plev)
  return(as.vector(t(outer(matnam,plev,paste,sep=':'))))
}

nxlevels <- function(n,f) {
  x <- attr(f,'x',exact=TRUE)
  plev <- rep(n,nlevels(f))
  if(is.null(x) || !is.matrix(x)) return(plev)
  nam <- attr(f,'xnam')
  matnam <- colnames(x)
  if(is.null(matnam)) matnam <- paste(nam,1:ncol(x),sep='') else matnam <- paste(nam,matnam,sep='')
  plev <- sub('.*:','',plev)
  return(as.vector(t(outer(matnam,plev,paste,sep=':'))))
}




