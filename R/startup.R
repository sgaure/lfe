# $Id: startup.R 1943 2016-04-07 23:08:38Z sgaure $

# set initial options, with defaults. Don't set if already set.
setoption <- function(...) {
  arg <- list(...)
  # don't change an existing option
  arg <- arg[sapply(paste('lfe',names(arg),sep='.'), function(n) is.null(getOption(n)))]
  # fill in from environment variable
  if(length(arg) == 0) return()
  nm <- names(arg)
  arg <- lapply(nm, function(n) {
    e <- Sys.getenv(paste('LFE',toupper(n),sep='_'))
    if(e != '') {
      val <- try(eval(parse(text=e)))
      if(inherits(val, 'try-error')) val <- arg[[n]]
      val
    } else arg[[n]]
  })
  names(arg) <- paste('lfe',nm,sep='.')
  do.call(options, arg)
}

.onLoad <- function(libname,pkgname) {
  setoption(usecg=FALSE, eps=1e-8, pint=1800L, accel=1L, bootmem=500, etol=c(1e-2,1e-12))

  if(is.null(cr <- getOption('lfe.threads'))) {
    cr <- as.integer(Sys.getenv('LFE_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_THREAD_LIMIT'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
    if(is.na(cr)) cr <- parallel::detectCores(all.tests=TRUE)
    if(is.na(cr) || cr < 1) {
      cr <- 1
    }
    options(lfe.threads=cr)
  }
}

.onUnload <- function(libpath) {
  options(lfe.usecg=NULL, lfe.eps=NULL,lfe.pint=NULL,lfe.accel=NULL,
          lfe.bootmem=NULL,lfe.threads=NULL, lfe.etol=NULL)
  library.dynam.unload('lfe',libpath)
}

if(!exists('anyNA')) anyNA <- function(x) any(is.na(x))

# Phase out cBind/rBind. Make sure we still work with older versions.
# Used in utils.R, chmethod.R and in some Rd-files (and vignettes)
# Make functions with (...) rather than just point to the right functions.
# The latter will raise a warning in checks about calls to .Internal functions.
RV <- R.Version()
rv <- paste(RV$major,RV$minor, sep='.')
mv <- sessionInfo()$otherPkgs$Matrix$Version
if(compareVersion('3.2-0',rv) > 0 || compareVersion('1.2-0',mv) > 0) {
  ..cbind.. <- quote(Matrix::cBind)
} else {
  ..cbind.. <- quote(cbind)
}
# make it local so we can avoid having the constant ..cbind.. in the namespace
mycbind <- local(function(...) {
  cl <- match.call()
  cl[[1L]] <- ..cbind..
  eval.parent(cl)
}, list(..cbind..=..cbind..))
rm(rv,RV,mv,..cbind..)


## numcores <- function() {
##   "This function is snatched from package 'multicore'"
##   systems <- list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
##                   freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
##                   linux = "grep processor /proc/cpuinfo 2>/dev/null|wc -l",
##                   irix = c("hinv |grep Processors|sed 's: .*::'", "hinv|grep '^Processor '|wc -l"),
##                   solaris = "/usr/sbin/psrinfo -v|grep 'Status of.*processor'|wc -l")
##   for (i in seq(systems)) if (length(grep(paste("^",
##                                                 names(systems)[i], sep = ""), R.version$os)))
##     for (cmd in systems[i]) {
##       a <- gsub("^ +", "", system(cmd, TRUE)[1])
##       if (length(grep("^[1-9]", a)))
##         return(as.integer(a))
##     }
##   NA
## }
