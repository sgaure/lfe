# function for doing multivariate non-linear tests on parameters
# tests the probability that fun(beta) > 0, where beta are the estimated parameters
# from felm()


#' Compute expectation of a function of the coefficients.
#' 
#' Use integration of the joint distribution of the coefficients to compute the
#' expectation of some function of the coefficients.  Can be used for
#' non-linear inference tests.
#' 
#' The function \code{nlexpect} integrates the function \code{fun(x)} over the
#' multivariate normal distribution specified by the point estimates and the
#' covariance matrix \code{vcov(est)}.  This is the expectation of
#' \code{fun(beta)} if we were to bootstrap the data (e.g. by drawing the
#' residuals anew) and do repeated estimations.
#' 
#' The list of coefficients used by \code{fun} must be specified in
#' \code{coefs}.
#' 
#' If the function is simple, it can be specified as a quoted expression like
#' \code{quote(a*b+log(abs(d)))}. In this case, if \code{coefs} is not
#' specified, it will be set to the list of all the variables occuring in the
#' expression which are also names of coefficients.
#' 
#' \code{fun} may return a vector of values, in which case a vector of
#' expectations is computed, like \code{quote(c(a*b, a^3-b))}. However, if the
#' expressions contain different variables, like \code{quote(c(a*b, d*e))}, a
#' quite compute intensive 4-dimensional integral will be computed, compared to
#' two cheap 2-dimensional integrals if you do them separately. There is nothing to gain
#' from using vector-valued functions compared to multiple calls to \code{nlexpect()}.
#' 
#' You may of course also integrate inequalites like \code{quote(abs(x1-0.2) >
#' 0.2)} to simulate the probability from t-tests or Wald-tests. See the
#' examples.
#' 
#' The function you provide will get an argument \code{...} if it does not have
#' one already.  It will also be passed an argument \code{.z} which contains
#' the actual coefficients in normalized coordinates, i.e. if \code{ch} is the
#' Cholesky decomposition of the covariance matrix, and \code{pt} are the point
#' estimates, the coefficients will be \code{pt + ch \%*\% .z}. The first argument
#' is a vector with names corresponding to the coefficients.
#'
#' If you specify \code{vectorized=TRUE}, your function will be passed a list with vectors
#' in its first argument. The function must
#' be able to handle a list, and must return a vector of the same length as the vectors
#' in the list.  If you pass an expression like \code{x < y}, each variable will be a vector.
#' If your function is vector valued, it must return a matrix where each
#' column is the values.
#' 
#' The \code{tol} argument specifies both the relative tolerance and the
#' absolute tolerance. If these should not be the same, specify \code{tol} as a
#' vector of length 2. The first value is the relative tolerance, the second is
#' the absolute tolerance. Termination occurs when at least one of the
#' tolerances is met.
#' 
#' The \code{...} can be used for passing other arguments to the integration
#' routine \code{cubature::cubintegrate} and the function to be integrated.
#' 
#' @param est object of class \code{"felm"} or \code{"lm"}, a result of a call to
#' \code{\link{felm}} or \code{lm}.
#' @param fun function of coefficients to be integrated. Can also be a
#' \code{quote}d expression.
#' @param coefs character. Names of coefficients to test. Only needed if
#' \code{fun} is a function.
#' @param ... other arguments passed to fun or the integration routine.
#' @param tol numeric. Tolerance for the computed integral.
#' @param lhs character. Name of the left hand side, if \code{est} has more
#' than one.
#' @param cv Covariance matrix to use in place of \code{vcov(est)}
#' @param istats logical. Should convergence information from the integration
#' routine be included as attributes?
#' @param flags list. Additional flags for the underlying integration routine. Not used after the
#' package \pkg{R2Cuba} disappeared.
#' @param max.eval integer. Maximum number of integral evaluations.
#' @param method character. A method specification usable by \code{cubature::cubintegrate}. 
#' The documentation there says that \code{"pcubature"} is good for smooth integrands of low dimensions.
#' @param vectorize logical or numeric. Use vectorized function evaluation from package
#' \pkg{cubature}. This can speed up integration significantly. If method is from the Cuba library
#' (i.e. not pcubature or hcubature), \code{vectorize} should be specified as a numeric, a vectorization
#' factor. The default is 128.
#' @return The function \code{nlexpect} computes and returns the expectation of
#' the function \code{fun(beta)}, with \code{beta} a vector of coefficients.
#' I.e., if the coefficients \code{beta} are bootstrapped a large number of
#' times, \code{nlexpect(est, fun)} should be equal to \code{mean(fun(beta))}.
#' @note An alternative to this method is to use the \code{bootexpr} argument
#' with \code{\link{felm}}, to do a Monte Carlo integration.
#' 
#' @seealso \code{\link{waldtest}}
#' @examples
#' 
#' N <- 100
#' x1 <- rnorm(N)
#' # make some correlation
#' x2 <- 0.1*rnorm(N) + 0.1*x1
#' y <- 0.1*x1 + x2 + rnorm(N)
#' summary(est <- felm(y ~ x1 + x2))
#' pt1 <- coef(est)['x1']
#' pt2 <- coef(est)['x2']
#' # expected values of coefficients, should match the summary
#' # and variance, i.e. square of standard errors in the summary
#' nlexpect(est, quote(c(x1=x1,x2=x2,var=c((x1-pt1)^2,(x2-pt2)^2))))
#' \donttest{
#' # the covariance matrix:
#' nlexpect(est, tcrossprod(as.matrix(c(x1-pt1,x2-pt2))))
#' }
#' #Wald test of single variable
#' waldtest(est, ~x1)['p.F']
#' # the same with nlexpect, i.e. probability for observing abs(x1)>abs(pt1) conditional
#' # on E(x1) = 0.
#' nlexpect(est, (x1-pt1)^2 > pt1^2, tol=1e-7, vectorize=TRUE)
#' # which is the same as
#' 2*nlexpect(est, x1*sign(pt1) < 0)
#'
#' # Here's a multivalued, vectorized example
#' nlexpect(est, rbind(a=x1*x2 < pt1, b=x1*x2 > 0), vectorize=TRUE)
#' \donttest{
#' 
#' # Non-linear test:
#'
#' # A simple one, what's the probability that product x1*x2 is between 0 and |E(x1)|?
#' nlexpect(est, x1*x2 > 0 & x1*x2 < abs(pt1), vectorize=TRUE)
#' # Then a more complicated one with the expected value of a polynomal in the coefficients
#' f <- function(x) c(poly=x[['x1']]*(6*x[['x1']]-x[['x2']]^2))
#' # This is the linearized test:
#' waldtest(est, f)['p.F']
#' # In general, for a function f, the non-linear Wald test is something like
#' # the following:
#' # expected value of function
#' Ef <- nlexpect(est, f, coefs=c('x1','x2'))
#' # point value of function
#' Pf <- f(c(pt1,pt2))
#' # similar to a Wald test, but non-linear:
#' nlexpect(est, function(x) (f(x)-Ef)^2 > Pf^2, c('x1','x2'), vectorize=TRUE)
#' # one-sided
#' nlexpect(est, function(x) f(x)-Ef > abs(Pf), c('x1','x2'), vectorize=TRUE)
#' # other sided
#' nlexpect(est, function(x) f(x)-Ef < -abs(Pf), c('x1','x2'), vectorize=TRUE)
#' }
#' 
#' @export nlexpect
nlexpect <- function(est, fun, coefs, ..., tol=getOption('lfe.etol'), lhs=NULL,
                     cv, istats=FALSE, flags=list(verbose=0), max.eval=200000L,
                     method=c('hcubature','pcubature','cuhre','suave','vegas','divonne'),
                     vectorize=FALSE) {

  mc <- match.call(expand.dots=FALSE)
  xargs <- names(mc[['...']])
  method <- match.arg(method)
  if(!requireNamespace('cubature', quietly=TRUE)) {
    warning('Package "cubature" not found.')
    return(NULL)
  }
#  adapt <- TRUE
#  adaptig <- switch(method,hcubature=cubature::hcubature,pcubature=cubature::pcubature,
#                    cuhre=cubature::cuhre,suave=cubature::suave, vegas=cubature::vegas)

  if(isTRUE(est$nostats) && missing(cv))
      stop('This test requires that felm() is run without nostats=TRUE; or specify a cv argument')


  # Find the covariance matrix
  if(missing(cv)) cv <- vcov(est, lhs=lhs)

  # Some kludge to be able to use non-quoted expression for fun
  afun <- substitute(fun)
  if(is.call(afun) || (is.name(afun) && (as.character(afun) %in% colnames(cv)))) {
    lfun <- as.list(afun)
    if(identical(lfun[[1]], quote(expression)))
        fun <- as.call(lfun[[2]])
    else if(!identical(lfun[[1]],quote(quote)) && !identical(lfun[[1]],quote(`function`)))
        fun <- afun
  }

  if(is.call(fun) || is.name(fun)) {
    # it's an expression. Figure out the coefficients used
    if(missing(coefs)) coefs <- intersect(all.vars(fun),colnames(cv))
    # make it a function
    fun <- local(function(x, ...) eval(fun,c(as.list(x),list(...))), list(fun=fun))
  } else if(is.function(fun)) {
    fa <- formals(fun)
#    nomatch <- !(xargs %in% names(fa)[-1])
#    if(any(nomatch))
#      warning('arguments ',paste(xargs[nomatch],collapse='/'), 
#              ' not among arguments in function to integrate')

    #add a ... formal if it doesn't exist
    if(!('.z' %in% names(fa)))
        formals(fun) <- c(fa,alist(.z=))
    if(!('...' %in% names(fa)))
        formals(fun) <- c(fa,alist(...=))
  }

  if(missing(coefs) || length(coefs %in% colnames(cv)) == 0)
      stop('No coefficients specified')
  # Find the coefficients
  cf <- drop(coef(est, lhs=lhs))[coefs]
  # and the Cholesky
  ch <- chol(cv[coefs,coefs,drop=FALSE])
  tch <- t(ch)
  # Now, we need to integrate fun(x) > 0 over the joint distribution of the parameters
  # We do this as follows. We integrate over a standard hypercube (-pi/2,pi/2) x (-pi/2,pi/2) x ...
  # adaptIntegrate can't take infinite limits.
  # We first transform these to (-Inf, Inf) with
  # z = tan(x)
  # the Jacobian determinant becomes the product of 
  # 1/cos(x)^2
  # We transform the integration variables with the covariance matrix to feed fun(),
  # then integrate fun(x) > 0 with the multivariate normal distribution.
  # we use the package cubature for the integration.
  K <- length(cf)
  if(is.numeric(vectorize)) nvec = if(vectorize<=1) 1 else vectorize
  if(isTRUE(vectorize)) nvec = 128 else nvec = 1
  if(nvec <= 1) {
    integrand <- function(x, ...) {
      jac <- prod(1/cos(x))^2
      z <- tan(x)
      # z is the standard normal (t really) multivariate
      #    dens <- prod(dnorm(z))
      dens <- prod(dt(z,est$df))
      beta <- drop(cf + tch %*% z)
      names(beta) <- coefs
      ret <- fun(beta, .z=z, ...)*jac*dens
      if(anyNA(ret)) stop('Function value is NA for argument: ',sprintf('%.2e ',beta))
      ret
    }
    sv <- fun(cf, .z=rep(0,K), ...)
  } else {
    integrand <- function(x, ...) {
      z <- tan(x)
      jac <- apply(x,2,function(x) prod(1/cos(x)))^2
      dens <- apply(z,2,function(zz) prod(dt(zz,est$df)))
      beta <- tch %*% z + cf
      lbeta <- vector('list',K)
      for(i in 1:nrow(beta)) lbeta[[i]] <- beta[i,]
      names(lbeta) <- coefs
      ret <- fun(lbeta, .z=z, ...)*jac*dens
      if(anyNA(ret)) stop('Function value is NA for argument: ',sprintf('%.2e ',beta))
      if(is.matrix(ret) && ncol(ret)==ncol(x)) ret else matrix(ret,ncol=ncol(x))
    }
    sv <- fun(as.list(cf), .z=matrix(0,K), ...)
  }
  

  fdim <- length(sv)
  if(length(tol) == 2) {
    reltol <- tol[1]
    abstol <- tol[2]
  } else {
    reltol <- abstol <- tol
  }

  eps <- 10*.Machine$double.eps
  ret <- cubature::cubintegrate(integrand,rep(-pi/2-eps,K),rep(pi/2-eps,K),fDim=fdim,
                      method=method,relTol=reltol,
                      absTol=abstol,maxEval=max.eval,nVec=nvec,...)
#  ret <- adaptig(integrand,rep(-pi/2-eps,K),rep(pi/2-eps,K),...,tol=reltol,
#                 absError=abstol,fDim=fdim,maxEval=max.eval,vectorInterface=vectorize)
  names(ret$integral) <- names(sv)
  if(is.array(sv)) {
    dim(ret$integral) <- dim(sv)
    dimnames(ret$integral) <- dimnames(sv)
  }
  if(!istats) return(ret$integral)
  
  names(ret)[match('integral',names(ret))] <- names(as.list(args(structure)))[1]  
  return(do.call(structure,ret))

#  ret <- R2ig(K,fdim,integrand,lower=rep(-1,K),...,
#              upper=rep(1,K),rel.tol=reltol,abs.tol=abstol,
#              flags=flags, max.eval=max.eval)
#  names(ret$value) <- names(sv)
#  if(is.array(sv)) {
#    dim(ret$value) <- dim(sv)
#    dimnames(ret$value) <- dimnames(sv)
#  }
#  if(ret$ifail != 0) warning('integration failed with: "',ret$message, '", use istats=TRUE to see details')
#  if(!istats) return(ret$value)
#  names(ret)[match('value',names(ret))] <- names(as.list(args(structure)))[1]
#  do.call(structure,ret)
}
