# return a data-frame with the group fixed effects, including zeros for references

#' Retrieve the group fixed effects
#' 
#' Compute the group fixed effects, i.e. the dummy parameters, which were swept
#' out during an estimation with \code{\link{felm}}.
#' 
#' For the case with two factors (the terms in the second part of the formula
#' supplied to \code{\link{felm}}), one reference in each connected component
#' is adequate when interpreting the results.
#' 
#' For three or more factors, no such easy method is known; for the
#' \code{"cholesky"} method- reference levels are found by analyzing the
#' pivoted Cholesky-decomposition of a slightly perturbed system.  The
#' \code{"kaczmarz"} method provides no rank-deficiency analysis, it is assumed
#' that the factors beyond the two first contribute nothing to the
#' rank-deficiency, so one reference in each is used.
#' 
#' If there are more than two factors, only the first two will be used to
#' report connected components.  In this case, it is not known which graph
#' theoretic concept may be used to analyze the rank-deficiency.
#' 
#' The standard errors returned by the Kaczmarz-method are bootstrapped,
#' keeping the other coefficients (from \code{\link{felm}}) constant, i.e. they
#' are from the variance when resampling the residuals.  If \code{robust=TRUE},
#' heteroskedastic robust standard errors are estimated. If \code{robust=FALSE}
#' and \code{cluster=TRUE}, clustered standard errors with the cluster
#' specified to \code{felm()} are estimated. If \code{cluster} is a factor, it
#' is used for the cluster definition.
#' 
#' @param obj object of class \code{"felm"}, usually, a result of a call to
#' \code{\link{felm}}
#' @param references a vector of strings.  If there are more than two factors
#' and you have prior knowledge of what the reference levels should be like
#' \code{references='id.23'}.  Not used with \code{method='kaczmarz'}
#' @param se logical.  Set to TRUE if standard errors for the group effects are
#' wanted.  This is \strong{very} time-consuming for large problems, so leave
#' it as FALSE unless absolutely needed.
#' @param method character string.  Either 'cholesky', 'cg', or the default
#' 'kaczmarz'.  The latter is often very fast and consumes little memory, it
#' requires an estimable function to be specified, see \code{\link{efactory}}.
#' The 'cholesky' method is no longer maintained as the author sees no use for
#' it.
#' @param ef function. A function of two variables, a vector of group fixed
#' effects and a logical, i.e. \code{function(v,addnames)}.  This function
#' should be estimable and is used to transform the raw-coefficients \code{v}
#' from the kaczmarz-method.  The second variable indicates whether the
#' function must return a named vector (if this is FALSE, one may skip the
#' names, saving memory allocations and time).
#' 
#' If a string is specified, it is fed to the \code{\link{efactory}}-function.
#' The default function is one which picks one reference in each component.
#' 
#' Can be set to \code{ef="ln"} to yield the minimal-norm solution from the
#' kaczmarz-method.
#' 
#' It can also be set to \code{ef="zm"} to get zero means (and intercept) in
#' one of the factors, and a reference in the other.
#' @param bN integer.  The number of bootstrap runs when standard errors are
#' requested.
#' @param robust logical. Should heteroskedastic standard errors be estimated?
#' @param cluster logical or factor. Estimate clustered standard errors.
#' @param lhs character vector. Specify which left hand side if \code{obj} has
#' multiple lhs.
#' @return The function \code{getfe} computes and returns a data frame
#' containing the group fixed effects.  It has the columns
#' \code{c('effect','se','obs','comp','fe','idx')}
#' 
#' \itemize{ \item \code{effect} is the estimated effect.  \item \code{se} is
#' the standard error.  \item \code{obs} is the number of observations of this
#' level.  \item \code{comp} is the graph-theoretic component number, useful
#' for interpreting the effects.  \item \code{fe} is the name of factor.  \item
#' \code{idx} is the level of the factor. }
#' 
#' With the Kaczmarz-method it's possible to specify a different estimable
#' function.
#' @keywords regression models
#' @examples
#' 
#' oldopts <- options(lfe.threads=2)
#' ## create covariates
#' x <- rnorm(4000)
#' x2 <- rnorm(length(x))
#' 
#' ## create individual and firm
#' id <- factor(sample(500,length(x),replace=TRUE))
#' firm <- factor(sample(300,length(x),replace=TRUE))
#' 
#' ## effects
#' id.eff <- rlnorm(nlevels(id))
#' firm.eff <- rexp(nlevels(firm))
#' 
#' ## left hand side
#' y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))
#' 
#' ## estimate and print result
#' est <- felm(y ~ x+x2 | id + firm)
#' summary(est)
#' ## extract the group effects
#' alpha <- getfe(est,se=TRUE)
#' 
#' ## find some estimable functions, with standard errors, we don't get
#' ## names so we must precompute some numerical indices in ef
#' idx <- match(c('id.5','id.6','firm.11','firm.12'),rownames(alpha))
#' alpha[idx,]
#' ef <- function(v,addnames) {
#'   w <- c(v[idx[[2]]]-v[idx[[1]]],v[idx[[4]]]+v[idx[[1]]],
#'          v[idx[[4]]]-v[idx[[3]]])
#'   if(addnames) names(w) <-c('id6-id5','f12+id5','f12-f11')
#'   w
#' }
#' getfe(est,ef=ef,se=TRUE)
#' options(oldopts)
#' \dontrun{
#' summary(lm(y ~ x+x2+id+firm-1))
#' }
#' 
#' @export getfe
getfe <- function(obj,references=NULL,se=FALSE,method='kaczmarz',ef='ref',bN=100, robust=FALSE, cluster=obj[['clustervar']], lhs=NULL) {

  if(length(obj$fe) == 0) return(NULL)
  if(!is.null(obj$numctrl) && obj$numctrl > 0)
      stop("Can't retrieve fixed effects when estimating with control variables")
  if(method == 'kaczmarz' || method == 'cg') {
    if(!is.null(references))
       warning('use estimable function (ef) instead of references in the Kaczmarz method')
    if(is.null(ef)) ef <- 'ln'
    if(!is.character(ef) && !is.function(ef))
      stop('ef must be a function when using the Kaczmarz method')
    if(method == 'cg') {
      oldopt <- options(lfe.usecg=TRUE)
      on.exit(options(oldopt))
    }
    return(getfe.kaczmarz(obj,se,ef=ef,bN=bN, robust=robust, cluster=cluster, lhs=lhs))
  }
  if(method != 'cholesky') stop('method must be either kaczmarz, cg, or cholesky')

  .Deprecated('',msg="Cholesky method is deprecated. Please consider using either 'kaczmarz' or 'cg'")
  attr(se,'sefactor') <- obj$sefactor
  attr(obj$fe,'references') <- references
  R <- obj$r.residuals
  # then the remaining.  This is usually sufficient.
  # we could also partition differently, just do the 'comp' adjustment accordingly
  # components
  ddlist <- makeddlist(obj$fe)
  gc()
  orignm <- attr(ddlist,'nm')
  comp <- 1
  res <- data.frame()
  for(dd in ddlist) {
#  res <- foreach(dd=ddlist,.combine=rbind,.init=data.frame()) %dopar% {
    dummies <- attr(dd,'dummies')
    keep <- attr(dd,'keep')
    comp <- attr(dd,'comp')
#    cat(date(),'comp dd',comp,'size',length(keep),'\n')
    Rhs <- as.vector(dummies %*% R[keep])
    names(Rhs) <- colnames(dd)
    alpha <- findfe(dd,Rhs,se)
    alpha[,'comp'] <- comp
    res <- rbind(res,alpha)
#    alpha
  }
  res <- res[orignm,]
  res[,'comp'] <- factor(res[,'comp'])    
# now, add factors telling which fe-group we're in
# the rownames are of the form <fe>.<idx>
  fefact <- strsplit(rownames(res),'.',fixed=TRUE)
  res[,'fe'] <- factor(unlist(lapply(fefact,function(l) l[[1]])))
  res[,'idx'] <- factor(unlist(lapply(fefact,function(l) paste(l[-1],collapse='.'))))
  return(res)
}
