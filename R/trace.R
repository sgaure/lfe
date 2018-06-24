# $Id: trace.R 1943 2016-04-07 23:08:38Z sgaure $
# compute the trace of a matrix.
# If we have a list of factors defining a projection, or a function for
# multiplying the matrix with a matrix, we use an iterative method
# from quantum physics, based on the general formula for the expectation
# of a quadratic form E(x' M x) = tr(MV) + E(x)' M E(x) where V=Cov(x). It reduces
# to 
# E(x' M x) = tr(M)
# when x has zero expectation and identity covariance matrix
# the most efficient (lowest variation) is to draw x uniformly from {-1,1}
# "Random phase vector for calculating the trace of a large matrix",
# T. Iitaka and T. Ebisuzaki, Physical Review E 69 (2004).







#' Compute trace of a large matrix by sample means
#' 
#' 
#' Some matrices are too large to be represented as a matrix, even as a sparse
#' matrix.  Nevertheless, it can be possible to compute the matrix vector
#' product fairly easy, and this is utilized to estimate the trace of the
#' matrix.
#' 
#' \code{mctrace} is used internally by \code{\link{fevcov}} and
#' \code{\link{bccorr}}, but has been made public since it might be useful for
#' other tasks as well.
#' 
#' For any matrix \eqn{A}, the trace equals the sum of the diagonal elements,
#' or the sum of the eigenvalues. However, if the size of the matrix is very
#' large, we may not have a matrix representation, so the diagonal is not
#' immediately available.  In that case we can use the formula \eqn{tr(A) =
#' E(x^t A x)}{tr(A) = E(x'Ax)} where \eqn{x} is a random vector with zero
#' expectation and \eqn{Var(x) = I}. We estimate the expecation with sample
#' means.  \code{mctrace} draws \eqn{x} in \eqn{\{-1,1\}^N}{{-1,1}}, and
#' evaluates \code{mat} on these vectors.
#' 
#' If \code{mat} is a function, it must be able to take a matrix of column
#' vectors as input.  Since \eqn{x^t A x = (Ax,x)}{x'Ax = (Ax,x)} is evaluated,
#' where \eqn{(\cdot,\cdot)}{(,)} is the Euclidean inner product, the function
#' \code{mat} can perform this inner product itself. In that case the function
#' should have an attribute \code{attr(mat,'IP') <- TRUE} to signal this.
#' 
#' If \code{mat} is a list of factors, the matrix for which to estimate the
#' trace, is the projection matrix which projects out the factors. I.e.  how
#' many dimensions are left when the factors have been projected out.  Thus, it
#' is possible to estimate the degrees of freedom in an OLS where factors are
#' projected out.
#' 
#' The tolerance \code{tol} is a relative tolerance.  The iteration terminates
#' when the normalized standard deviation of the sample mean (s.d. divided by
#' absolute value of the current sample mean) goes below \code{tol}.  Specify a
#' negative \code{tol} to use the absolute standard deviation.  The tolerance
#' can also change during the iterations; you can specify
#' \code{tol=function(curest) {...}} and return a tolerance based on the
#' current estimate of the trace (i.e. the current sample mean).
#' 
#' @param mat square matrix, Matrix, function or list of factors.
#' @param N integer. if \code{mat} is a function, the size of the matrix is
#' specified here.
#' @param tol numeric. Tolerance.
#' @param maxsamples numeric. Maximum number of samples in the expectation
#' estimation.
#' @param trname character. Arbitrary name used in progress reports.
#' @param init numeric. Initial guess for the trace.
#' @return An estimate of the trace of the matrix represented by \code{mat} is
#' returned.
#' @examples
#' 
#'   A <- matrix(rnorm(25),5)
#'   fun <- function(x) A %*% x
#'   sum(diag(A))
#'   sum(eigen(A,only.values=TRUE)$values)
#'   # mctrace is not really useful for small problems.
#'   mctrace(fun,ncol(A),tol=0.05)
#'   # try a larger problem (3000x3000):
#'   f1 <- factor(sample(1500,3000,replace=TRUE))
#'   f2 <- factor(sample(1500,3000,replace=TRUE))
#'   fl <- list(f1,f2)
#'   mctrace(fl,tol=-5)
#'   # exact:
#'   length(f1) - nlevels(f1) - nlevels(f2) + nlevels(compfactor(fl))
#' 
#' @export mctrace
mctrace <- function(mat, N, tol=1e-3, maxsamples=Inf,
                    trname='', init) {
  if(is.matrix(mat) || inherits(mat,'Matrix')) {
    return(structure(sum(diag(mat)), sd=0, iterations=0))
  } else if(is.list(mat) && all(sapply(mat, is.factor))) {
    N <- length(mat[[1]])
    fun <- function(v,trtol) colSums(demeanlist(v, mat)*v)
  } else if(!is.function(mat)) {
    stop('mat must be function, factor list or matrix')
  } else {
    if(missing(N)) stop('N (vector length) must be specified with mat a function')
    if(isTRUE(attr(mat,'IP'))) {
      # inner product is done by the function itself.
      fun <- mat
    } else {
      fun <- function(v,trtol) colSums(mat(v)*v)
    }
  }

  if(!is.function(tol)) eps <- function(x) tol else eps <- tol
  
  threads <- getOption('lfe.threads')
  if(maxsamples < threads) threads <- maxsamples
  maxB <- getOption("lfe.bootmem") * 1e+06
  maxvpt <- maxB %/% (2*8*N*threads)
  if(maxvpt*threads > 4096) maxvpt <- 4096 %/% threads
#  if(maxvpt == 0) {
#    maxvpt <- 1
#  }
  # ensure at least 8 vectors in first iteration
  if(threads >= 8)
      vpt <- 1
  else
      vpt <- 16 %/% (threads+1)

  
  blk <- vpt*threads
  if(blk > maxsamples) blk <- maxsamples
  i <- 0
  tr <- 0
  sqsum <- 0
  NN <- 0
  last <- start <- Sys.time()
  # get a clue about the tolerance.
#  cureps <- eps(as.numeric(fun(0, trtol=0)))/2
  if(missing(init)) init <- N
  cureps <- eps(init)/2
  pint <- getOption('lfe.pint')

  while(NN < maxsamples && (NN < 8 ||
                            (cureps > 0 && relsd > cureps) ||
                            (cureps < 0 && sd > -cureps))) {
    i <- i+1
    now <- Sys.time()
    if(NN > 0) {
      remaining <- as.integer((Ntarget-NN)/(NN/as.numeric(now-start, units='secs')))
      if(remaining > pint && as.numeric(now - last, units='secs') > pint) {
        message('  *** trace ',trname,' sample ',NN,' of ',Ntarget,
                ',mean ',signif(tr/NN,3), ', sd ',signif(sd,3),', target ', signif(cureps,3),
                ', expected finish at ',
                now + remaining) 
        last <- now
      }
    }
    
    ests <- fun(matrix(sample(c(-1,1), N*blk, replace=TRUE), N), trtol=abs(cureps))
    gc()
#    cat('ests : ', mean(ests), ' ', sd(ests)); print(fivenum(ests))
    NN <- NN + blk
    tr <- tr + sum(ests)
    sqsum <- sqsum + sum(ests^2)
    rm(ests)
    # compute sd for the mean tr/NN. It's sqrt(E(x^2) - E(X)^2)/sqrt(NN)
    if(NN > 1)
        sd <- sqrt(sqsum/(NN-1) - tr^2/((NN*(NN-1))))/sqrt(NN)
    else
        sd <- 0

#    message(trname,' sd=',sd,' relsd=',relsd,' NN=',NN, ' cureps=',cureps)
    cureps <- eps(tr/NN)

    # try to figure out how many iterations are needed to obtain the
    # desired tolerance.
    sdtarget <- if(cureps < 0) -cureps else cureps*abs(tr/NN)
    if(NN == 1) sd <- sqrt(2)*sdtarget
    relsd <- sd/abs(tr/NN)
    Ntarget <- as.integer((sd/sdtarget)^2*NN)
    if(is.na(Ntarget)) stop('Too much variance in trace samples for ',trname, ', sd ',sd,', cureps ',cureps)
    vpt <- 1 + (Ntarget-NN) %/% threads
    if(vpt > maxvpt) vpt <- maxvpt
    if(vpt < 1) vpt <- 1
    blk <- vpt*threads
  }
#  if(last > start) cat('\n')
  dt <- as.numeric(Sys.time()-start, units='secs')
  if(dt > pint)
      message('  *** trace ',trname,' ',NN, ' samples finished in ', as.integer(dt), ' seconds') 
  structure(tr/NN, sd=sd, iterations=NN)
}
