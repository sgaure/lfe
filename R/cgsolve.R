# $Id: cgsolve.R 1990 2016-04-14 13:48:36Z sgaure $
# From "A practical termination criterion for the conjugate gradient method", E.F. Kaasschieter
# BIT 28 (1988), 308-322  (2.6). 
# return the vector phi_i(x) for i=1:j

phi <- function(j, x, alpha, beta) {
  if(j == 1) return(1)

  psiprev <- 1
  phiprev <- 1-alpha[1]*x
  res <- c(1,phiprev)
  if(j == 2) return(res)

  for(i in 3:j) {
    psi <- phiprev + beta[i-1]*psiprev
    phi <- phiprev - alpha[i]*x*psi
    res <- c(res,phi)
    phiprev <- phi
    psiprev <- psi
  }
  as.numeric(res)
}

newmu <- function(oldmu, i, alpha, beta) {
  if(all(phi(i,oldmu, alpha,beta) > 0)) return(oldmu)
  u <- 1e-3
  y <- 0
  z <- oldmu
  # binary search
  x <- oldmu
  while(z - y > u*y) {
    x <- (y+z)/2
    if(all(phi(i,x,alpha,beta) > 0))
        y <- x
    else
        z <- x
  }
  return(x)
}

newmus <- function(oldmus, i, alpha, beta) {
  res <- NULL
  for(j in seq_along(oldmus)) {
    res <- c(res, newmu(oldmus[j],i,alpha[,j],beta[,j]))
  }
  res
}

#stop('debug')
# conjugate gradient solver Ax = b, b may be matrix
# If A is a function, it must be able to take a matrix argument
# Algorithm 3 from Kaasschieter (1988)






#' Solve a symmetric linear system with the conjugate gradient method
#' 
#' 
#' \code{cgsolve} uses a conjugate gradient algorithm to solve the linear
#' system \eqn{A x = b} where \eqn{A} is a symmetric matrix.  \code{cgsolve} is
#' used internally in \pkg{lfe} in the routines \code{\link{fevcov}} and
#' \code{\link{bccorr}}, but has been made public because it might be useful
#' for other purposes as well.
#' 
#' 
#' The argument \code{A} can be a symmetric matrix or a symmetric sparse matrix
#' inheriting from \code{"Matrix"} of the package \pkg{Matrix}.  It can also be
#' a function which performs the matrix-vector product. If so, the function
#' must be able to take as input a matrix of column vectors.
#' 
#' If the matrix \code{A} is rank deficient, some solution is returned.  If
#' there is no solution, a vector is returned which may or may not be close to
#' a solution.  If \code{symmtest} is \code{FALSE}, no check is performed that
#' \code{A} is symmetric. If not symmetric, \code{cgsolve} is likely to raise
#' an error about divergence.
#' 
#' The tolerance \code{eps} is a relative tolerance, i.e.  \eqn{||x - x_0|| <
#' \epsilon ||x_0||} where \eqn{x_0} is the true solution and \eqn{x} is the
#' solution returned by \code{cgsolve}. Use a negative \code{eps} for absolute
#' tolerance.  The termination criterion for \code{cgsolve} is the one from
#' \cite{Kaasschieter (1988)}, Algorithm 3.
#' 
#' Preconditioning is currently not supported.
#' 
#' If \code{A} is a function, the test for symmetry is performed by drawing two
#' random vectors \code{x,y}, and testing whether \eqn{|(Ax, y) - (x, Ay)| <
#' 10^{-6} sqrt((||Ax||^2 + ||Ay||^2)/N)}, where \eqn{N} is the vector length.
#' Thus, the test is neither deterministic nor perfect.
#' 
#' @param A matrix, Matrix or function.
#' @param b vector or matrix of columns vectors.
#' @param eps numeric. Tolerance.
#' @param init numeric. Initial guess.
#' @param symmtest logical. Should the matrix be tested for symmetry?
#' @param name character. Arbitrary name used in progress reports.
#' @return
#' 
#' A solution \eqn{x} of the linear system \eqn{A x = b} is returned.
#' @seealso \code{\link{kaczmarz}}
#' @references Kaasschieter, E. (1988) \cite{A practical termination criterion
#' for the conjugate gradient method}, BIT Numerical Mathematics,
#' 28(2):308-322.  \url{http://link.springer.com/article/10.1007\%2FBF01934094}
#' @examples
#' 
#'   N <- 100000
#' # create some factors
#'   f1 <- factor(sample(34000,N,replace=TRUE))
#'   f2 <- factor(sample(25000,N,replace=TRUE))
#' # a matrix of dummies, which probably is rank deficient
#'   B <- makeDmatrix(list(f1,f2))
#'   dim(B)
#' # create a right hand side
#'   b <- as.matrix(B %*% rnorm(ncol(B)))
#' # solve B' B x = B' b
#'   sol <- cgsolve(crossprod(B), crossprod(B, b), eps=-1e-2)
#'   #verify solution
#'   sqrt(sum((B %*% sol - b)^2))
#' 
#' @export cgsolve
cgsolve <- function(A, b, eps=1e-3,init=NULL, symmtest=FALSE, name='') {
  start <- Sys.time()
  precon <- attr(A,'precon')
  oneini <- 0
  if(is.matrix(b) || inherits(b,'Matrix')) N <- nrow(b) else N <- length(b)
  if(is.matrix(A) || inherits(A,'Matrix') ) {
    if(symmtest && !isSymmetric(A)) stop("matrix is not symmetric")
    fun <- function(v) as.matrix(A %*% v) 
  } else if(is.function(A)) {
    fun <- function(v) as.matrix(A(v))
    if(symmtest) {
      # A symmetric matrix satisfies (Ax,y) = (Ay,x) for every y and x
      x <- as.matrix(runif(N,0.7,1.4))
      y <- as.matrix(runif(N,0.7,1.4))
      Ax <- fun(x)
      Ay <- fun(y)
      if(abs(sum(Ax * y) - sum(Ay * x)) > 1e-6*sqrt(mean(Ax^2) + mean(Ay^2)))
          warning(deparse(substitute(A)), ' seems to be non-symmetric')
      rm(x,y,Ax,Ay)
    }
  } else {
    stop("A must be a matrix, Matrix or a function")
  }

  b <- as.matrix(b)
  start <- last <- Sys.time()

#  if(is.null(init)) init <- matrix(rnorm(nrow(b)*ncol(b)),nrow(b))
  if(is.null(init)) {
    x <- matrix(0,N,ncol(b))
    r <- b
  } else {
    x <- as.matrix(init)
    if(ncol(x) == 1 && ncol(b) > 1)
        x <- matrix(init, N, ncol(b))
    r <- b - fun(x)
  }
  p <- r
  tim <- 0
# first iteration outside loop
  oldr2 <- colSums(r^2)
  minr2 <- sqrt(max(oldr2))
  Ap <- fun(p)
  alpha <- oldr2/colSums(p * Ap)
  mu <- 1/alpha
  allalpha <- matrix(alpha,1)

  x <- x + t(alpha*t(p)) 
  r <- r - t(alpha*t(Ap))
  r2 <- colSums(r^2)
  
  res <- matrix(0,nrow(x),ncol(x))
  origcol <- 1:ncol(b)
  k <- 0
  kk <- 0
  allbeta <- NULL
  if(length(eps) == 1) eps <- rep(eps,ncol(b))
  negeps <- eps < 0
  eps <- abs(eps)
  eps <- ifelse(negeps, eps, eps/(1+eps))
  pint <- getOption('lfe.pint')
  while(TRUE) {

    k <- k+1
    delta <- mu * ifelse(negeps,1,sqrt(colSums(x^2)))

    now <- Sys.time()
    dt <- as.numeric(now-last, units='secs')
    if(dt > pint) {
      last <- now
      message(date(), ' CG iter ',k,' ',name,' target=',signif(min(eps),3),
              ' delta=',signif(max(sqrt(r2)/delta),3), ' vecs=',length(delta), ' mu ',signif(min(mu),3))
    }


    done <- sqrt(r2) < eps * delta | (k >= N)

    if(any(done)) {
#      message('finished vec ',ilpretty(origcol[done]))
      # remove the finished vectors
      res[,origcol[done]] <- x[,done,drop=FALSE]
      origcol <- origcol[!done]
      eps <- eps[!done]
      negeps <- negeps[!done]
      x <- x[,!done, drop=FALSE]
      if(length(origcol) == 0) break;
      p <- p[,!done, drop=FALSE]
      r <- r[,!done, drop=FALSE]
      r2 <- r2[!done]
      oldr2 <- oldr2[!done]
      mu <- mu[!done]
      allalpha <- allalpha[,!done, drop=FALSE]
      allbeta <- allbeta[,!done, drop=FALSE]
      kk <- 0
    }
    r2rms <- sqrt(max(r2))
    minr2 <- min(minr2,r2rms)
    if(( kk > 1000 && r2rms > 100*minr2) || (k > 100 && r2rms > 10000*minr2)) {
      warning('cgsolve (',name,') seems to diverge, iter=',k,', ||r2||=',r2rms,
              ' returning imprecise solutions')
      res[,origcol[seq_len(ncol(x))]] <- x
      return(res)
    }

    kk <- kk + 1
    beta <- r2/oldr2
    allbeta <- rbind(allbeta, beta)
# we have created some C-functions to handle some operations.
# Solely to avoid copying large matrices.
    p <- .Call(C_pdaxpy,r,p,beta)
#    p <- r + t(beta*t(p))
    Ap <- fun(p)
#    gc()
    alpha <- r2/.Call(C_piproduct, Ap, p)
#    alpha <- r2/colSums(Ap*p)
#    message(name, ' ', k, ' alpha '); print(alpha)
    allalpha <- rbind(allalpha, alpha)
    x <- .Call(C_pdaxpy, x, p, alpha)
    r <- .Call(C_pdaxpy, r, Ap, -alpha)
    oldr2 <- r2
    r2 <- colSums(r^2)
    mu <- newmus(mu, k+1, allalpha, allbeta)
  }
#  message('CG iters:',k)
  now <- Sys.time()
  dt <- as.numeric(now-start, units='secs')
  if(dt > pint)
      message('  *** cgsolve(',name,') finished with ',k,' iters in ',as.integer(dt),' seconds')
  res
}
