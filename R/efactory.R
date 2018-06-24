# $Id: efactory.R 1943 2016-04-07 23:08:38Z sgaure $






#' Create estimable function
#' 
#' Creates an estimable function for a factor-structure.
#' 
#' There are several possibilities for the input parameter \code{opt}.
#' \itemize{ \item \code{"ref"} yields an estimable function which is similar
#' to the default one in \code{\link{lm}}, one reference is forced to \code{0}
#' in each connected component.  \item \code{"zm"} Similar to \code{"ref"}, but
#' the factor which does not contain a reference is made to have zero mean, and
#' an intercept is added.  \item \code{"zm2"} Similar to \code{"zm"}, but both
#' factors are made to have zero mean.  \item \code{"ln"} Least norm function.
#' This will yield the raw coefficients from the Kaczmarz-method, i.e. the
#' solution with smallest norm. This function is not estimable.  } Note that in
#' the case with more than two factors, it is not known how to analyze the
#' factors to find the structure of the rank-deficiencies, i.e. the estimable
#' functions.  In this case, the factors beyond the first two are assumed not
#' to contribute to the rank-deficiency beyond a single dimension in each.
#' Both \code{"ref"} and \code{"zm"} keep one such reference at zero in each of
#' these factors.  This is the common method when using dummies.
#' 
#' In the case that interactions are specified in the model, i.e. with
#' \code{x:f} in the second part of the formula, these terms are not analyzed
#' to create an estimable function. Only the pure \code{f} terms are used for
#' this purpose.  It is assumed that the \code{x:f} terms are all identified.
#' Note that in this case, all the levels of \code{f} are included.
#' 
#' @param obj object of class \code{"felm"}, usually, a result of a call to
#' \code{\link{felm}}.
#' @param opt character.  Which type of estimable function.
#' @param ... various.
#' @return A function of two parameters \code{function(v,addnames)}.  An
#' estimable function (i.e. the result is the vector of some length \code{N})
#' of the input vector \code{v}. When \code{addnames==TRUE} the returned vector
#' should have names, and optionally an attribute \code{"extra"} which is a
#' list of vectors of length \code{N} which may be used to code additional
#' information.
#' @note The author is open to suggestions for other estimable functions, i.e.
#' other useful normalizations of the solutions.
#' 
#' It is not strictly necessary that the \code{obj} argument is of class
#' \code{"felm"}, any list with entries \code{"fe"} and \code{"cfactor"} of the
#' appropriate form will do. That is, \code{list(fe=fl,cfactor=compfactor(fl))}
#' where \code{fl} is the list of factors defining the component structure.
#' I.e. if the model is \code{y ~ ... |id + firm}, we have
#' \code{fl=list(id=id,firm=firm)}.
#' @examples
#' 
#' oldopts <- options(lfe.threads=1)
#' id <- factor(sample(5000,50000,replace=TRUE))
#' firm <- factor(sample(3000,50000,replace=TRUE))
#' fl <- list(id=id,firm=firm)
#' obj <- list(fe=fl,cfactor=compfactor(fl))
#' ## the trivial least-norm  transformtion, which by the way is non-estimable
#' print(ef <- efactory(obj,'ln'))
#' is.estimable(ef,fl)
#' ## then the default
#' print(ef <- efactory(obj,'ref'))
#' is.estimable(ef,fl)
#' # get the names of the coefficients, i.e. the nm-variable in the function
#' head(evalq(nm,environment(ef)))
#' options(oldopts)
#' 
#' @export efactory
efactory <- function(obj, opt='ref', ...) {

  # only factors without covariates are relevant to analyze
  purefes <- sapply(obj$fe, function(f) is.null(attr(f,'x',exact=TRUE)))
  pfe <- obj$fe[purefes]

#  allnm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='\003')))
  allnm <- unlist(lapply(names(obj$fe),function(n) xlevels(n,obj$fe[[n]],sep='\003')))
  
# the names of the dummies, e.g. id.4 firm.23
#  nm <- unlist(lapply(names(pfe),function(n) paste(n,levels(pfe[[n]]),sep='\003')))
  nm <- unlist(lapply(names(pfe),function(n) xlevels(n,pfe[[n]],sep='\003')))
  # create an index where the pure fe's belong in the full array
  allpos <- match(nm,allnm)

  mkallvec <- function(x) {res <- rep(NA,length(allnm)); res[allpos] <- allpos[x]; res;}
# how many obervations for each level
  lobs <- lapply(pfe,table)
  obs <- unlist(lobs)  
  names(obs) <- unlist(lapply(names(lobs), function(n) paste(n,'\003',names(lobs[[n]]),sep='')))

#  allobs <- unlist(lapply(obj$fe,table))
  allobs <- unlist(lapply(obj$fe,function(f) {
    x <- attr(f,'x',exact=TRUE)
    if(is.null(x)) return(table(f))
    if(!is.matrix(x)) return(table(f))
    return(rep(table(f), ncol(x)))
  }))  

  if(length(pfe) == 2) {
    # now, find the component of each parameter, i.e. each level.  We do this
    # by finding the first occurence of each level, i.e. match(levels(f),f)
    comp <- factor(unlist(lapply(pfe, function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
  } else if(length(pfe) > 2) {
    # we should formally assign unique component numbers for factors beyond the second
    comp <- factor(unlist(lapply(pfe[1:2], function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
    exlvls <- (nlevels(comp)+1):(nlevels(comp)+1 + length(pfe)-3)
    comp <- as.factor(c(comp,unlist(mapply(rep,exlvls,unlist(lapply(pfe[3:length(pfe)],nlevels))))))
  } else {
    comp <- factor(rep(1,length(obs)))
    ncomp <- 1
  }

  if(length(pfe) == 0) {
    nm <- unlist(lapply(names(obj$fe),function(n) xlevels(n,obj$fe[[n]],sep='.')))
    return(function(v,addnames) {
      if(!addnames) return(v)
      names(v) <- nm
      v
    })
  }

  refnames <- unlist(tapply(obs,comp,function(l) names(which.max(l))))
  # now v[refnames] will be the reference values

  refno <- match(refnames,nm)
  refsub <- refno[comp]
  # refsub is a vector, in entry i it contains the reference for entry i
  # i.e. we should do a v <-  v - v[refsub]
  # subtract all references.
  # but for the main components we should only do this for the
  # factor in which the references is.
  # for the other factor we should add the reference
  # thus we need two versions of refsub, one with NA's in the
  # reference factor, one with NA's in the other, then we must
  # replace NA's with zero before subtracting
  # so which ones belong to which factor?
  # make a factor to decide

  fef <- factor(unlist(lapply(names(pfe),function(n) rep(n,nlevels(pfe[[n]])))))
#  allfef <- factor(unlist(lapply(names(obj$fe),function(n) rep(n,nlevels(obj$fe[[n]])))))
  allfef <- factor(unlist(lapply(names(obj$fe),function(n) nxlevels(n,obj$fe[[n]]))))

  # level of the factor
#  idx <- factor(unlist(lapply(obj$fe,function(f) levels(f))))
  idx <- factor(unlist(lapply(obj$fe,function(f) {
    x <- attr(f,'x',exact=TRUE)
    if(is.null(x) || !is.matrix(x)) return(levels(f))
    return(rep(levels(f), ncol(x)))
  })))
  # then figure out in which factor the reference is
  # make sure to allow '.' in factor names
  rf <- sub('(^.*)\003..*$','\\1',refnames)
  # now, create a refsubs which is the ones to be subtracted
  # each refsub belonging to somthing else than the reference factor
  # should be NA'ed.
  if(length(pfe) > 2) {
    extra <- (length(refno)-length(pfe)+3):length(refno)
    sw <- c(names(pfe)[c(2,1)],rep('.NA',length(pfe)-2))
  } else {
    swap <- if(length(pfe) == 2) c(2,1) else 1
    sw <- names(pfe)[swap]
    extra <- integer(0)
  }
  names(sw) <- names(pfe)
  otherf <- sw[rf]
  # which should we not subtract?  
  # Those which are 
  nosub <- fef != rf[comp]
  refsubs <- refsub
  refsubs[nosub] <- NA

  # which should we add, those which are different from the reference factor
  noadd <- fef != otherf[comp]
  refsuba <- refsub
  refsuba[noadd] <- NA

  extrarefs <- refno[extra]

  # now, what if we should centre on the means?
  # there are two variants, either centre on the means in
  # both factors, or only in the one without a reference
  # we create a factor zfact which describes the groups
  # to mean.  

  
  # create a minimal environment for a function

  extrarefs <- allpos[extrarefs]
  refsubs <- mkallvec(refsubs)
  refsuba <- mkallvec(refsuba)
  obs <- allobs
  fef <- allfef
#  nm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='.')))
  nm <- unlist(lapply(names(obj$fe),function(n) xlevels(n,obj$fe[[n]],sep='.')))
#  allcomp <- rep(0,sum(sapply(obj$fe,nlevels)))
  allcomp <- rep(0,length(allnm))
  allcomp[allpos] <- comp
  comp <- allcomp
  fenv <- list(extrarefs=extrarefs,refsubs=refsubs,refsuba=refsuba,fef=fef,nm=nm)


  ef <- switch(as.character(opt),
         ln={
           local(function(v,addnames) {
             if(addnames) {
               names(v) <- nm
               attr(v,'extra') <- list(obs=obs,comp=comp,fe=fef,idx=idx)
             }
             v
           },list(obs=obs,comp=comp,fe=fef,idx=idx,nm=nm))
         },
         # one reference in each component
         ref={
           fenv$comp <- comp
           fenv$mkallvec <- mkallvec
           local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             if(addnames) {
               names(v) <- nm
               attr(v,'extra') <- list(obs=obs,comp=comp,fe=fef,idx=idx)
             }
             v
           },fenv)
         },
         zm={
           # now, what if we want zero-means on the other-factor?
           # we will then get an intercept for each component, and
           # zero means.  It's those which are in nosub, but partitioned
           # into components. We may do this faster now that we've
           # separated it from the ordinary 'ref'
           zfact <- comp
           zfact[!nosub] <- NA

           enames <- paste('icpt',1:ncomp,sep='.')
           zcomp <- factor(c(comp,1:ncomp))
           oo <- order(zcomp)
           fenv$oo <- oo
           fenv$zfact <- zfact
           fenv$zcomp <- zcomp[oo]
           fenv$enames <- enames
           fenv$obs <- c(obs,table(obj$cfactor))[oo]
           ef <- local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             means <- tapply(v,zfact,mean)
             mn <- means[zfact]
             mn <- ifelse(is.na(mn),0,mn)
             v <- v - mn
             v <- c(v,means)[oo]
             if(addnames) {
               names(v) <- c(nm,enames)[oo]
               attr(v,'extra') <- list(obs=obs,comp=zcomp)
             }
             v
           },fenv)
         },
         # one reference in each component, but zero-means in the other factor
         # and an intercept
         zm2={
           # both factors (but not the extra factors):
           # the interaction between comp and fef forms these groups,
           zfact <- interaction(comp,fef)
           # but skip the extra factors
           zfact[as.integer(comp) > ncomp] <- NA
           zfact <- factor(zfact)
           # and the means should be added to the intercepts
           # i.e. from the same component, add two and two
           # ifact should consist of the components of the
           # levels of zfact.
           ifact <- factor(as.integer(gsub('^([0-9]+).*','\\1',levels(zfact))),exclude=NA)
           enames <- paste('icpt',1:ncomp,sep='.')
           zcomp <- factor(c(comp,1:ncomp))
           oo <- order(zcomp)
           fenv$oo <- oo
           fenv$zfact <- zfact
           fenv$zcomp <- zcomp[oo]
           fenv$enames <- enames
           fenv$obs <- c(obs,table(obj$cfactor))[oo]
           fenv$ifact <- ifact
           ef <- local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             means <- tapply(v,zfact,mean)
             mn <- means[zfact]
             mn <- ifelse(is.na(mn),0,mn)
             v <- v - mn
             icpt <- tapply(means,ifact,sum)
             v <- c(v,icpt)[oo]
             if(addnames) {
               names(v) <- c(nm,enames)[oo]
               attr(v,'extra') <- list(obs=obs,comp=zcomp)
             }
             v
           },fenv)
         },

         stop(paste('estimable function',opt,'not recognized'))
         )

# try to byte compile the stuff
  ef <- compiler::cmpfun(ef,list(optimize=3))
  if(length(pfe) <= 2 && as.character(opt) != 'ln' && all(purefes)) 
    attr(ef,'verified') <- TRUE
  ef
} 
