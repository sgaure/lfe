#' Centre vectors on multiple groups
#'
#' @description
#'  Uses the method of alternating projections to centre
#'  a (model) matrix on multiple groups, as specified by a list of factors.
#'  This function is called by \code{\link{felm}}, but it has been
#'  made available as standalone in case it's needed. In particular, if
#' one does not need transformations provided by R-formulas but have the covariates present
#' as a matrix or a data.frame, a substantial amount of time can be saved in the centering.
#' @export
#\usage{
#demeanlist(mtx, fl, icpt=0, eps=getOption('lfe.eps'),
#           threads=getOption('lfe.threads'),
#           progress=getOption('lfe.pint'),
#           accel=getOption('lfe.accel'),
#           randfact=TRUE,
#           means=FALSE,
#           weights=NULL,
#           scale=TRUE,
#           .inplace=FALSE)
#}

#' @param mtx matrix whose columns form vectors to be group-centred. mtx
#'  can also be a list of vectors or matrices, such as a data frame.
#' @param fl list of factors defining the grouping structure.
#' @param icpt the position of the intercept, this column is removed from
#' the result matrix.
#' @param eps a tolerance for the centering.
#' @param threads an integer specifying the number of threads to use.
#' @param progress integer. If positive, make progress reports (whenever a
#' vector is centered, but not more often than every \code{progress} minutes).
#' @param accel integer. Set to 1 if Gearhart-Koshy acceleration should be done.
#' @param randfact logical. Should the order of the factors be randomized?
#' This may improve convergence.
#' @param means logical. Should the means instead of the demeaned matrix be
#'  returned? Setting \code{means=TRUE} will return \code{mtx -  demeanlist(mtx,...)},
#' but without the extra copy.
#' @param weights numeric. For weighted demeaning.
#' @param scale logical. Specify scaling for weighted demeaning.
#' @param na.rm logical which indicates what should happen when the data 
#' contain \code{NA}s. If TRUE, rows in the input \code{mtx} are removed
#' prior to centering. If FALSE, they are kept, leading to entire groups becoming NA
#' in the output. 
#' @param attrs list. List of attributes which should be attached to the output. 
#' Used internally.
#'
#' @details
#' For each column \code{y} in \code{mtx}, the equivalent of the
#' following centering is performed, with \code{cy} as the result.
#' \preformatted{  
#' cy <- y; oldy <- y-1
#' while(sqrt(sum((cy-oldy)**2)) >= eps) {
#'  oldy <- cy
#'  for(f in fl) cy <- cy - ave(cy,f)
#' }
#' }
#'
#' Each factor in \code{fl} may contain an
#' attribute \code{'x'} which is a numeric vector of the same length as
#' the factor. The centering is then not done on the means of each group,
#' but on the projection onto the covariate in each group.  That is, with a
#' covariate \code{x} and a factor \code{f}, it is like projecting out the
#' interaction \code{x:f}.  The \code{'x'} attribute can also be a matrix of column
#' vectors, in this case it can be beneficial to orthogonalize the columns,
#' either with a stabilized Gram-Schmidt method, or with the simple
#' method \code{x \%*\% solve(chol(crossprod(x)))}.
#'
#' The \code{weights} argument is used if a weighted projection is
#' computed.  If \eqn{W} is the diagonal matrix with \code{weights} on the
#' diagonal, \code{demeanlist} computes \eqn{W^{-1} M_{WD} W x} where \eqn{x} is
#' the input vector, \eqn{D} is the matrix of dummies from \code{fl} and
#' \eqn{M_{WD}} is the projection on the orthogonal complement of
#' the range of the matrix \eqn{WD}.  It is possible to implement the
#' weighted projection with the \code{'x'} attribute mentioned above, but
#' it is a separate argument for convenience.
#' If \code{scale=FALSE}, \code{demeanlist} computes \eqn{M_{WD} x} without
#' any \eqn{W} scaling.
#' If \code{length(scale) > 1}, then \code{scale[1]} specifies whether
#' the input should be scaled by \eqn{W}, and \code{scale[2]} specifies
#' whether the output should be scaled by \eqn{W^{-1}}.  This is just
#' a convenience to save some memory copies in other functions in the package.
#' 
#' Note that for certain large datasets the overhead in \code{\link{felm}}
#' is large compared to the time spent in \code{demeanlist}. If the data
#' are present directly without having to use the formula-interface to
#' \code{felm} for transformations etc, it is possible to run
#' \code{demeanlist} directly on a matrix or \code{"data.frame"} and do the
#' OLS "manually", e.g. with something like
#' \code{cx <- demeanlist(x,...); beta <- solve(crossprod(cx), crossprod(cx,y))}
#' 
#' In some applications it is known that a single centering iteration is
#' sufficient. In particular, if \code{length(fl)==1} and there is no
#' interaction attribute \code{x}.  In this case the centering algorithm is
#' terminated after the first iteration. There may be other cases, e.g. if
#' there is a single factor with a matrix \code{x} with orthogonal columns. If
#' you have such prior knowledge, it is possible to force termination after
#' the first iteration by adding an attribute \code{attr(fl, 'oneiter') <-
#' TRUE}.  Convergence will be reached in the second iteration anyway, but
#' you save one iteration, i.e. you double the speed.
#' 
#' @return
#' If \code{mtx} is a matrix, a matrix of the same shape, possibly with 
#' column \code{icpt} deleted.
#' 
#' If \code{mtx} is a list of vectors and matrices, a list of the same
#' length is returned, with the same vector and matrix-pattern, but the
#' matrices have the column \code{icpt} deleted.  
#' 
#' If \code{mtx} is a \code{'data.frame'}, a \code{'data.frame'} 
#' with the same names is returned; the \code{icpt} argument is ignored.
#' 
#' If \code{na.rm} is specified, the return value has an attribute \code{'na.rm'} with a vector of
#' row numbers which has been removed. In case the input is a matrix or list, the same rows
#' are removed from all of them. Note that removing NAs incurs a copy of the input, so if
#' memory usage is an issue and many runs are done, one might consider removing NAs from the data set entirely.
#' @note
#' The \code{accel} argument enables Gearhart-Koshy acceleration as
#' described in Theorem 3.16 by Bauschke, Deutsch, Hundal and Park in "Accelerating the
#' convergence of the method of alternating projections",
#' Trans. Amer. Math. Soc. 355 pp 3433-3461 (2003).
#'
#' \code{demeanlist} will use an in place transform to save memory, provided the \code{mtx}
#' argument is unnamed. Thus, as always in R, you shouldn't use temporary variables
#' like \code{tmp <- fun(x[v,]); bar <- demeanlist(tmp,...); rm(tmp)}, it will be much better to
#' do \code{bar <- demeanlist(fun(x[v,]),...)}. However, demeanlist allows a construction like
#' \code{bar <- demeanlist(unnamed(tmp),...)} which will use an in place transformation, i.e. tmp
#' will be modified, quite contrary to the usual semantics of R.
#' 
#' @examples
#' oldopts <- options(lfe.threads=1)
#' ## create a matrix
#' mtx <- data.frame(matrix(rnorm(999),ncol=3))
#' # a list of factors
#' rgb <- c('red','green','blue')
#' fl <- replicate(4, factor(sample(rgb,nrow(mtx),replace=TRUE)), simplify=FALSE)
#' names(fl) <- paste('g',seq_along(fl),sep='')
#' # centre on all means
#' mtx0 <- demeanlist(mtx,fl)
#' head(data.frame(mtx0,fl))
#' # verify that the group means for the columns are zero
#' lapply(fl, function(f) apply(mtx0,2,tapply,f,mean))
#' options(oldopts)
demeanlist <- function(mtx,fl,icpt=0L,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads'),
		       progress=getOption('lfe.pint'),
                       accel=getOption('lfe.accel'),
                       randfact=TRUE,
                       means=FALSE,
                       weights=NULL,
                       scale=TRUE,
                       na.rm=FALSE,
                       attrs=NULL) {

  if(length(fl) == 0) {
    if(means) {
      foo <- unlist(utils::as.relistable(mtx))
      foo[] <- 0
      return(utils::relist(foo))
    }
    return(eval.parent(substitute(mtx)))
#    return(mtx)
  }

  # Here we used to just .Call(C_demeanlist, mtx, fl,...)
  # However, we have rewritten this part because C_demeanlist checks the NAMED
  # status of mtx, to avoid copying if possible. But mtx is certainly NAMED > 0 inside
  # this function.  So instead, we evaluate C_demeanlist in the parent environment, without
  # touching the actual arguments in here
  mf <- match.call()
  ff <- formals(sys.function())
  # This is the argument sent to C_demeanlist, in this order:
  ff <- ff[match(c("mtx","fl","icpt","eps","threads","progress","accel","means","weights","scale","attrs"),
             names(ff), 0L)]
  m <- names(ff)[names(ff) %in% names(mf)]
  ff[m] <- mf[m]
  # make a new environment to store our reordered fl, enclosed by our caller's frame
  # This is just to avoid clutter in case of error messages
  env <- new.env(parent=parent.frame())
  assign('C_demeanlist',C_demeanlist,envir=env)
  assign('unnamed',unnamed,envir=env)
  assign('.fl',eval.parent(ff[['fl']]), envir=env)
  .fl <- NULL # avoid check warning
  ff[['fl']] <- quote(.fl)
  if(randfact && length(get('.fl',env)) > 2) {
    # This is delayed for the bizarre reason that before we rewrote this code
    # mtx was forced before reordering of fl, so we just keep it that way to avoid
    # changing the test output.
    delayedAssign('..fl',.fl[order(runif(length(.fl)))],env,env)
    ff[['fl']] <- quote(..fl)
  }

  # Then some NA-handling, at the request of David Hugh-Jones, June 11, 2018
  badrows <- NULL
  delist <- FALSE
  isDT <- FALSE
  if(isTRUE(na.rm)) {
    # we need to touch and copy the mtx and fl
    # mtx is either a vector, a matrix or a list (data.frame/data.table) of such
    mtx <- eval.parent(ff[['mtx']])
    if(!is.list(mtx)) {mtx <- list(mtx); delist <- TRUE}
    # find bad rows in any of the elements of mtx
    badrows <- NULL
    for(i in seq_along(mtx)) {
      object <- mtx[[i]]
      d <- dim(object)
      if (length(d) > 2L)  next
      bad <- seq_along(object)[is.na(object)]
      if (length(bad) == 0L) next
      if (length(d)) {
        bad <- unique(((bad - 1)%%d[1L]) + 1L)
      } 
      badrows <- union(badrows,bad)
    }
    if(length(badrows) > 0) {
      isDF <- is.data.frame(mtx)
      newmtx <- lapply(mtx,function(m) if(length(dim(m)) > 1) m[-badrows,,drop=FALSE] else m[-badrows])
      rm(mtx)
      if(isDF) newmtx <- as.data.frame(newmtx)
      if(delist) newmtx <- newmtx[[1]]
      fl <- eval(ff[['fl']],env)
      N <- length(fl[[1]])
      newfl <- lapply(fl, function(f) factor(f[-badrows]))
      assign('.mtx',newmtx,envir=env)
      assign('.fl',newfl,envir=env)
      ff[['fl']] <- quote(.fl)
      rm(newfl)
      ff[['attrs']] <- c(ff[['attrs']],list(na.rm=sort(badrows)))
    } else {
      assign('.mtx',mtx,envir=env)
    }
    ff[['mtx']] <- quote(unnamed(.mtx))
  }
  eval(as.call(c(list(quote(.Call), quote(C_demeanlist)), ff)), env)
}

