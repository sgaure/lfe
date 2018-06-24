# $Id: compfactor.R 1943 2016-04-07 23:08:38Z sgaure $






#' Find the connected components
#' 
#' 'compfactor' computes the connected components of the dummy-part of the
#' model.
#' 
#' If there are more than two factors and \code{WW=FALSE}, only the first two
#' will be used.
#' 
#' If \code{WW=TRUE} and \code{length(fl) > 2}, the component structure will be
#' as in "A Note on the Determination of Connectedness in an N-Way Cross
#' Classification" by D.L. Weeks and D.R. Williams, Technometrics, vol 6 no 3,
#' August 1964. I.e. in each component, the coefficients within each factor are
#' comparable, that is, their difference is estimable even if there are more
#' than two factors.  That is, one may use one reference in each factor in each
#' component, and interpret the coefficients within a component as usual. This
#' is not an exhaustion of all the estimable functions.  There is somewhat more
#' about this in one of the vignettes.
#' 
#' @param fl a list of factors defining the dummies
#' @param WW logical. Use Weeks and Williams components
#' @return A factor of the same length as the factors in the input argument.
#' It defines the connected components. E.g. \code{nlevels(compfactor(fl))}
#' will yield the number of connected components.
#' @examples
#' 
#' ## create two factors
#' f1 <- factor(sample(300,400,replace=TRUE))
#' f2 <- factor(sample(300,400,replace=TRUE))
#' 
#' ## find the components
#' cf <- compfactor(list(f1=f1,f2=f2))
#' 
#' ## show the third largest component
#' fr <- data.frame(f1,f2,cf)
#' fr[cf==3,]
#' 
#' @export compfactor
compfactor <- function(fl, WW=FALSE) {
  if(length(fl) == 0) return(factor(NULL))
  N <- length(fl[[1]])
  purefls <- sapply(fl,function(f) is.null(attr(f,'x',exact=TRUE)))
  fl <- fl[purefls]
  if(length(fl) <= 1) return(factor(rep(1,N)))
  if(WW && length(fl) > 2) {
    cf <- factor(.Call(C_wwcomp,fl))
  } else {
    cf <- factor(.Call(C_conncomp,fl[1:2]))
  }
  cf
}
