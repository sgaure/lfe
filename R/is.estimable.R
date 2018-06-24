
#' Verify estimability of function
#' 
#' Verify that a function you have written for \code{\link{getfe}} is indeed
#' estimable.
#' 
#' When writing custom estimable functions for \code{\link{getfe}}, the
#' function \code{is.estimable} can be used to test it for estimability.
#' \code{is.estimable()} solves the sparse residual system with the Kaczmarz
#' method, using two different initial values. Then \code{ef()} is applied to
#' the two solutions. If the value of \code{ef()} differs by more than
#' \code{1e-5} in any coordinate, \code{FALSE} is returned, otherwise
#' \code{TRUE} is returned.  If \code{keepdiff=TRUE}, the vector of differences
#' is attached as an attribute \code{'diff'} to the returned logical value.  If
#' you have problems with estimability, it is a fair guess that those entries
#' with a difference in absolute values smaller than, say, \code{1e-5} are
#' estimable, whereas the others are not.
#' 
#' @param ef function.  The function to be verified.
#' @param fe list of factors.
#' @param R numeric.  Vector of residuals, if \code{NULL}, a random one is
#' created.
#' @param nowarn logical. Set to \code{TRUE} if \code{is.estimable} should not
#' throw a warning for non-estimable functions.
#' @param keepdiff logical. Return differences between two different runs of
#' the Kaczmarz method.
#' @param threshold numeric. Threshold for determining estimability.
#' @return Returns a logical.
#' @seealso \code{\link{getfe}}
#' @examples
#' 
#' oldopts <- options(lfe.threads=1)
#' ## create individual and firm
#' id <- factor(sample(5000,50000,replace=TRUE))
#' firm <- factor(sample(3000,50000,replace=TRUE))
#' 
#' ## create some estimable functions. It's faster to
#' ## use numerical indices in ef rather than strings, and the input v
#' ## to ef has no names, we have to add them when requested
#' ef <- function(v,addnames) {
#'   w <- c(v[6]-v[5],v[7000]+v[5],v[7000]-v[6000])
#'   if(addnames) names(w) <-c('id6-id5','f2k+id5','f2k-f1k')
#'   w
#' }
#' is.estimable(ef,list(id=id,firm=firm))
#' 
#' ## Then make an error; in the last coordinate, sum two firms
#' ef <- function(v,addnames) {
#'   w <- c(v[6]-v[5],v[7000]+v[5],v[7000]+v[6000])
#'   if(addnames) names(w) <-c('id6-id5','f2k+id5','f2k-f1k')
#'   w
#' }
#' is.estimable(ef, list(id=id,firm=firm), keepdiff=TRUE)
#' options(oldopts)
#' 
#' @export is.estimable
is.estimable <- function(ef,fe,R=NULL,nowarn=FALSE,keepdiff=FALSE, threshold=500*getOption('lfe.eps')) {
  if(!is.function(ef)) stop('ef must be a function')
  N <- sum(unlist(lapply(fe,function(f) {
    x <- attr(f,'x', exact=TRUE)
    if(is.matrix(x)) nlevels(f)*ncol(x) else nlevels(f)
  })))

  if(is.null(R)) {
    # make a suitable residual

    nr <- length(fe[[1]])
    vec <- unlist(lapply(fe,function(f) {
      x <- attr(f,'x', exact=TRUE)
      if(is.matrix(x)) return(unlist(apply(x,2,function(cl) cl*runif(nlevels(f))[f])))
      r <- runif(nlevels(f))[f]
      if(is.null(x)) r else unlist(r*x)
    }))
    dim(vec) <- c(nr, length(vec)/nr)
    R <- rowSums(vec)
  }
  v1 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  v2 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  df <- max(abs(v1-v2))
  if(df > threshold) {
    bad <- which.max(abs(v1-v2))
    badname <- names(bad)
    if(!nowarn)
      warning('non-estimable function, largest error ',
              format(df,digits=1),' in coordinate ',bad, ' ("',badname,'")')
    return(structure(FALSE,diff=if(!keepdiff) NULL else v1-v2))
  }
  structure(TRUE,diff=if(!keepdiff) NULL else v1-v2)
}
