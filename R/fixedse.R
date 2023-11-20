#' Compute standard errors for fixed effects
#'
#' fixedse computes the standard errors for the fixed effects when there is only one.
#' While [getfe()] can provide standard errors, it does so by bootstrapping
#' for general estimable functions. In the special case that there's only a single fixed
#' effect, and the estimable function is just the levels, this function can be used to
#' compute the fixed effects without bootstrapping. It requires that [felm()]
#' is run with keepX=TRUE.
#' @name fixedse
#' @param est 'felm' object. The result of a call to [felm()].
#' @param E Matrix. Estimable function. Not used at the moment.
#' @param lhs character. Name of the left hand side, if more than one.
#' @return numeric. Vector of standard errors.
#' @examples
#' x <- rnorm(1000)
#' f <- factor(sample(5, 1000, replace = TRUE))
#' y <- x + (1:5)[f] + rnorm(1000)
#' est <- felm(y ~ x | f, keepX = TRUE)
#' # both bootstrap and computed se:
#' cbind(getfe(est, ef = efactory(est, "ref"), se = TRUE), fse = fixedse(est))
#' # compare with lm:
#' summary(lm(y ~ x + f - 1))
#' @export
#' @keywords internal
#' @importFrom Matrix colSums rowSums solve
fixedse <- function(est, lhs = NULL, E) {
  if (length(est$fe) == 0) {
    return(numeric())
  }
  if (!is.null(est$clustervar)) {
    stop("fixedse() does not work with clustered standard errors")
  }
  s2 <- sum(residuals(est, lhs = lhs)^2) / df.residual(est)

  if (length(est$fe) != 1) {
    stop("fixedse() only works for a single fixed effect")
    # Now, we must have ef an estimable function
    if (missing(E) || !inherits(E, "Matrix")) stop("Must specify estimable function with more than one factor")
    D <- makeDmatrix(est$fe)
    DE <- D %*% E
    dd <- crossprod(DE)
    dig <- diag(solve(dd))
    #    dig <- sampdiag(dd,eps,ncol(E))
    DX <- crossprod(DE, est$X)
    cvc <- chol(vcov(est, lhs = lhs))
    return(structure(sqrt(s2 * dig + colSums(tcrossprod(cvc, solve(dd, DX))^2)),
      names = colnames(E)
    ))
  }
  f <- est$fe[[1]]
  # Use the Woodbury identity on (D' M_X D)^{-1} = (D' (I-X(X'X)^{-1}X')) D)^{-1}
  # to obtain (D'D)^{-1}(D'D + D'X(X' M_D X)^{-1})X'D) (D'D)^{-1}
  # where (X' M_D X)^{-1} is the unscaled vcov(est)
  if (is.null(est$X) || ncol(est$X) == 0) {
    sqrt(s2 / as.vector(table(f)))
  } else {
    sqrt(s2 / as.vector(table(f)) + colSums(tcrossprod(
      chol(vcov(est, lhs = lhs)),
      crowsum(est$X, f, mean = TRUE)
    )^2))
  }
}

sampdiag <- function(A, eps = 0.01, K) {
  # Start with 10 vectors
  NN <- 0
  N <- 100
  if (is.matrix(A) || inherits(A, "Matrix")) {
    K <- ncol(A)
  } else if (missing(K)) {
    stop("Must specify K=ncol(A) when A is function")
  }

  meansd <- 2 * eps
  sums <- sqsums <- numeric(K)
  cgeps <- sqrt(K) * eps / 10
  while (meansd > eps) {
    message("meanSD = ", meansd, " ", NN, " new N is ", N)
    N <- as.integer(max(1, min(1024 * 1024 * getOption("lfe.bootmem") / (3 * K * 8), N)))
    #    N <- min(N,500L)
    vv <- matrix(sample(c(1, -1), N * K, replace = TRUE), K)
    vAv <- vv * solve(A, vv)
    #    vAv <- vv*cgsolve(A, vv, cgeps)
    rm(vv)
    sums <- sums + rowSums(vAv)
    sqsums <- sqsums + rowSums(vAv^2)
    rm(vAv)
    NN <- NN + N
    sd <- sqrt((sqsums - sums^2 / NN) / (NN - 1)) / sqrt(NN)
    res <- sums / NN
    relsd <- sd / (0 + abs(res))
    print(fivenum(relsd))
    meansd <- max(relsd) # max(relsd[relsd > eps])
    #    meansd <- max(abs(sd))
    if (is.na(meansd)) break
    # How many Ns do we need to reach eps?
    N <- NN * ((meansd / eps)^2 - 1)
  }
  message("NN is ", NN)
  res
}
