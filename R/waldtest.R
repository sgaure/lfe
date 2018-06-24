#' Compute Wald test for joint restrictions on coefficients
#' 
#' Compute a Wald test for a linear hypothesis on the coefficients.  Also
#' supports Delta-approximation for non-linear hypotheses.
#' 
#' The function \code{waldtest} computes a Wald test for the H0: R beta = r,
#' where beta is the estimated vector \code{coef(object)}.
#' 
#' If \code{R} is a character, integer, or logical vector it is assumed to
#' specify a matrix which merely picks out a subset of the coefficients for
#' joint testing. If \code{r} is not specified, it is assumed to be a zero
#' vector of the appropriate length.
#' 
#' \code{R} can also be a formula which is linear in the estimated
#' coefficients, e.g. of the type \code{~Q-2|x-2*z} which will test the joint
#' hypothesis Q=2 and x=2*z.
#' 
#' If \code{R} is a function (of the coefficients), an approximate Wald test
#' against H0: \code{R(beta) == 0}, using the Delta-method, is computed.
#' 
#' In case of an IV-estimation, the names for the endogenous variables in
#' \code{coef(object)} are of the type \code{"`Q(fit)`"} which is a bit dull to
#' type; if all the endogenous variables are to be tested they can be specified
#' as \code{"endovars"}. It is also possible to specify an endogenous variable
#' simply as \code{"Q"}, and \code{waldtest} will add the other syntactic sugar
#' to obtain \code{"`Q(fit)`"}.
#' 
#' The \code{type} argument works as follows. If \code{type=='default'} it is
#' assumed that the residuals are i.i.d., unless a cluster structure was
#' specified to \code{\link{felm}}. If \code{type=='robust'}, a heteroscedastic
#' structure is assumed, even if a cluster structure was specified in
#' \code{\link{felm}}.
#' 
#' @param object object of class \code{"felm"}, a result of a call to
#' \code{\link{felm}}.
#' @param R matrix, character, formula, function, integer or logical.
#' Specification of which exclusions to test.
#' @param r numerical vector.
#' @param type character. Error structure type.
#' @param lhs character. Name of left hand side if multiple left hand sides.
#' @param df1 integer. If you know better than the default df, specify it here.
#' @param df2 integer. If you know better than the default df, specify it here.
#' @return The function \code{waldtest} computes and returns a named numeric
#' vector containing the following elements.
#' 
#' \itemize{ \item \code{p} is the p-value for the Chi^2-test \item \code{chi2}
#' is the Chi^2-distributed statistic.  \item \code{df1} is the degrees of
#' freedom for the Chi^2 statistic.  \item \code{p.F} is the p-value for the F
#' statistics \item \code{F} is the F-distributed statistic.  \item \code{df2}
#' is the additional degrees of freedom for the F statistic. }
#' 
#' The return value has an attribute \code{'formula'} which encodes the
#' restrictions.
#' @seealso \code{\link{nlexpect}}
#' @examples
#' 
#' x <- rnorm(10000)
#' x2 <- rnorm(length(x))
#' y <- x - 0.2*x2 + rnorm(length(x))
#' #Also works for lm
#' summary(est <- lm(y ~ x + x2  ))
#' # We do not reject the true values
#' waldtest(est, ~ x-1|x2+0.2|`(Intercept)`)
#' # The Delta-method coincides when the function is linear:
#' waldtest(est, function(x) x - c(0, 1, -0.2))
#' 
#' @export waldtest
waldtest <- function(object, R, r, type=c('default','iid','robust','cluster'), lhs=NULL, df1, df2) {
  if(inherits(object,'felm') && object$nostats) stop('No Wald test for objects created with felm(nostats=TRUE)')

  # We make a chi^2 to test whether the equation R theta = r holds.
  # The chi^2 is computed according to Wooldridge (5.34, 10.59).
  # I.e. a Wald test W = N*(beta' (R V^{-1} R')^{-1} beta) where beta = R theta - r
  # W is chi2 with length(r) df,
  # and V is th covariance matrix.

  # First, find V. It's in either object$vcv, object$robustvcv or object$clustervcv
  if(is.null(lhs) && length(object$lhs) > 1) {
    stop('Please specify lhs=[one of ',paste(object$lhs, collapse=','),']')
  }
  if(!is.null(lhs) && is.na(match(lhs, object$lhs)))
      stop('Please specify lhs=[one of ',paste(object$lhs, collapse=','),']')

  type <- type[1]
  if(identical(type,'default')) {
    if(is.null(object$clustervar))
        V <- vcov(object, type='iid', lhs=lhs)
    else
        V <- vcov(object, type='cluster', lhs=lhs)
  } else
      V <- vcov(object, type=type, lhs=lhs)

#  if(is.null(lhs) && length(object$lhs) == 1) lhs <- object$lhs
  cf <- coef(object)
  if(is.matrix(cf))
      nmc <- rownames(cf)
  else
      nmc <- names(cf)

  if(inherits(R,'formula') || is.call(R) || is.name(R)) {
    Rr <- formtoR(R, nmc)
    R <- Rr[,-ncol(Rr), drop=FALSE]
    r <- Rr[,ncol(Rr)]
  } else if(is.function(R)) {
    # non-linear stuff. Compute value and gradient of R
  if(!requireNamespace('numDeriv', quietly=TRUE)) {warning("package numDeriv must be available to use non-linear Wald test"); return(NULL)}
    pt <- coef(object,lhs=lhs)
    pt[is.na(pt)] <- 0
    val <- R(pt)
    if(is.null(dim(val))) dim(val) <- c(length(val), 1)
    gr <- numDeriv::jacobian(R,pt)
    if(is.null(dim(gr))) dim(gr) <- c(1,length(gr))
  } else  if(!is.matrix(R)) {
    # it's not a matrix, so it's a list of parameters, either
    # names, logicals or indices
    if(is.null(R)) R <- nmc
    if(is.character(R)) {
      ev <- match('endovars', R)
      if(!is.na(ev)) {
        # replace with endogenous variables
        R <- c(R[-ev],object$endovars)
      }
      # did user specify any of the endogenous variables?
      fitvars <- paste('`',R,'(fit)`',sep='')
      fitpos <- match(fitvars,nmc)
      # replace those which are not NA
      noNA <- which(!is.na(fitpos))
      R[noNA] <- fitvars[noNA]
      Ri <- match(R, nmc)
      if(anyNA(Ri)) stop("Couldn't find variables ",paste(R[is.na(Ri)],collapse=','))
      R <- Ri
    } else if(is.logical(R)) {
      R <- which(R)
    }
    # here R is a list of positions of coefficients
    # make the projection matrix.

    RR <- matrix(0,length(R),length(coef(object,lhs=lhs)))
    for(i in seq_along(R)) {
      RR[i,R[i]] <- 1
    }
    R <- RR
  } 
  # Two cases here. If R is a function, we do a non-linear delta test against 0, otherwise
  # we do the ordinary Wald test
  if(is.function(R)) {
    W <- as.numeric(t(val) %*% solve(gr %*% V %*% t(gr)) %*% val)
    if(missing(df1)) df1 <- length(val)
  } else {
    if(missing(r) || is.null(r))
        r <- rep(0,nrow(R))
    else if(length(r) != nrow(R)) stop('nrow(R) != length(r)')
    cf <- coef(object, lhs=lhs)
    cf[is.na(cf)] <- 0
    beta <- R %*% cf - r
    V[is.na(V)] <- 0   # ignore NAs
    W <- try(sum(beta * solve(R %*% V %*% t(R),beta)), silent=TRUE)
    
    if(inherits(W,'try-error')) 
        W <- as.numeric(t(beta) %*% pinvx(R %*% V %*% t(R)) %*% beta)
  }
  # W follows a chi2(Q) distribution, but the F-test has another
  # df which is ordinarily object$df. However, if there are clusters
  # the df should be reduced to the number of clusters-1
  if(missing(df2)) {
    df2 <- object$df
    if((!is.null(object$clustervar) && type %in% c('default','cluster')) ) {
      df2 <- min(nlevels(object$clustervar[[1]])-1, df2)
    }
  }

  if(missing(df1))
      df1 <- length(beta)

  F <- W/df1
  # F follows a F(df1,df2) distribution
  if(is.function(R)) frm <- R else  frm <- Rtoform(R,r,nmc)
  
  structure(c(p=pchisq(W, df1, lower.tail=FALSE), chi2=W, df1=df1,
              p.F=pf(F,df1,df2, lower.tail=FALSE), F=F, df2=df2),
            formula=frm)
}

# convert a formula which is a set of linear combinations like ~x+x3 | x2-x4+3 to
# matrices R and r such that R %*% coefs = r
# the vector r is return as the last column of the result
formtoR <- function(formula, coefs) {

  conv <- function(f) formtoR(f, coefs)

  lf <- as.list(formula)
  if(lf[[1]] == as.name('~') || lf[[1]] == as.name('quote')) return(conv(lf[[2]]))
  # here we have a single formula w/o '~' in front, e.g. x+x3|x2-x4, or just x+x3
  # split off parts '|' in a loop
  R <- NULL
#  if(length(lf) != 1) stop('length of ',lf, ' is != 1')
#  lf <- as.list(lf[[1]])
  op <- lf[[1]]
  if(op == as.name('|')) {
    return(rbind(conv(lf[[2]]), conv(lf[[3]])))
  } else if(op == as.name('+')) {
    if(length(lf) == 2) return(conv(lf[[2]])) # unary +
    return(conv(lf[[2]]) + conv(lf[[3]]))
  } else if(op == as.name('-')) {
    if(length(lf) == 2) return(-conv(lf[[2]])) # unary -
    return(conv(lf[[2]]) - conv(lf[[3]]))    
  } else if(op == as.name('*')) {
    f1 <- conv(lf[[2]])
    f2 <- conv(lf[[3]])
    # the first one must be a numeric, i.e. only last column filled in
    # and it's negative
    fac <- -f1[length(f1)]
    return(fac * conv(lf[[3]])) 
  } else if(is.name(op)) {
    res <- matrix(0,1,length(coefs)+1)
    pos <- match(as.character(op), coefs)
    if(is.na(pos)) {
      ivspec <- paste("`",as.character(op),"(fit)`", sep='')
      pos <- match(ivspec, coefs)
    }
    if(is.na(pos)) stop("Can't find ", op, " among coefficients ", paste(coefs, collapse=','))
    res[pos] <- 1
    return(res)
  } else if(is.numeric(op)) {
    return(matrix(c(rep(0,length(coefs)), -op), 1))
  } else {
    stop('Unkwnown item ',as.character(op), ' in formula ',formula)
  }
}

Rtoform <- function(R,r, coefs) {
  coefs <- gsub('`','',coefs,fixed=TRUE)
  form <- paste('~',paste(apply(R, 1, function(row) {
    w <- which(row != 0)
    rw <- paste(' ', row[w], '*`', coefs[w], '`', collapse=' + ', sep='')
    rw <- gsub('+ -',' - ',rw, fixed=TRUE)
    rw <- gsub(' 1*','',rw, fixed=TRUE)
    rw <- gsub('(fit)','',rw, fixed=TRUE)
    rw
  }), ' + ', -r, collapse='|', sep=''))
  form <- gsub('+ -','-',form, fixed=TRUE)
  form <- gsub(' 0.',' .',form, fixed=TRUE)
  form <- gsub('+ 0','',form, fixed=TRUE)
  local(as.formula(form))
}
