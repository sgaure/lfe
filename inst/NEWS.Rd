\name{NEWS}
\title{lfe news}
\section{Changes in version 2.8-5}{
  \itemize{
    \item Corrected the vignette about estimable functions. The 
    change in the \code{base::sample()} function made the examples meaningless, and I
    forgot to fix it in 2.8-3.
  }
}
\section{Changes in version 2.8-4}{
  \itemize{
    \item Support for nested clusters and reghdfe-option by Grant McDermott and Karl
    Dunkle Werner.
    \item Changed some extern declarations in the C-code to comply with
    standards enforced by the upcoming gcc10.
    \item Changed inner workings of the \code{unnamed()} and
    \code{demeanlist()} functions so they
    still conserve some memory with the new reference count machinery.
  }
}

\section{Changes in version 2.8-3}{
  \itemize{
    \item Fixed a bug in Crowsum which prevented computation of clustered standard
    errors in large datasets (rows x covariates exceeds 2
    billion). Thanks to mcmalone for the report.
    \item Minor updates to documentation and some cosmetic code changes.
    \item Updates to configure.ac at the request of Brian Ripley
    \item Updates to tests to reflect the new uniform \code{sample()} algorithm
    \item summary.felm now defaults to robust errors if
    \code{getOption(lfe.robust=TRUE)}. As with the other options, lfe.robust is
    initialized from the environment variable LFE_ROBUST upon load.
    \item A felm-object now works as expected with \code{formula()}.
    \item A felm-object now works as expected with \code{sandwich::bread()}.
  }
}
\section{Changes in version 2.8}{
  \itemize{
    \item Added functionality for closed form standard errors on the fixed
    effects if standard normal and only one fixed effect.
  }
}
\section{Changes in version 2.7}{
  \itemize{
    \item Added a \code{na.rm} argument to \code{demeanlist} to remove
    rows with NAs in the input data prior to centering. Previously these were silently set to 0.

    \item A bug in the negative Eigenvalue adjustment for multiway
    clustering has been fixed.

    \item Dependence on package \pkg{R2Cuba} has been removed because
    that package is orphaned and will be removed from CRAN.
  }
}
\section{Changes in version 2.6}{
\itemize{
  \item If running with robust or clustered standard errors,
  \code{getfe} now adds a column 'se' which is a copy of the
  'clusteredse' or 'robustse'. This to ensure compatibility with package
  \pkg{broom}.

  \item Added a routine \code{sargan()} to compute Sargan's S.
}
}
\section{Changes in version 2.5}{
  \itemize{
  \item \code{demeanlist()} will now do the transformation in place if
  given an unnamed input. This saves memory in various other places.

  \item Multiway clustered standard errors now incorporate the negative
  Eigenvalue adjustment of Cameron, Gelbach, Miller. I.e. they are
  forced to zero. Can be switched off by \code{felm(...,psdef=FALSE)}.

  \item Fixed a bug in \code{btrap} which could yield wrong bootstrapped standard
  errors in \code{getfe}.

  \item A function \code{nlexpect()} has been added to integrate
  functions of estimated coefficients over the residual distribution.

  \item \code{waldtest()} now supports non-linear restrictions with the Delta-method.
  }
}

\section{Changes in version 2.4 (2015-10-26}{
\itemize{
  \item A bug in bias correction with clustered errors has been fixed.

  \item A bug in the handling of expressions f:g with two factors has
  been fixed.

  \item A function \code{estfun.felm} has been introduced to improve
  interoperability with the package \pkg{multiwayvcov}.  Also a function
  \code{weights.felm} has been introduced.

  \item A bug in \code{residuals.felm} has been fixed.

  \item \code{felm} now supports bootstrapping an arbitrary expression
  involving estimated coefficients. This can also be done by repeated
  calls to \code{felm}, but the builtin functionality avoids repeated
  calls to \code{model.frame} and \code{model.matrix} which can be quite
  time-consuming for large datasets.
}
}

\section{Changes in version 2.3 (2015-05-25)}{
  \itemize{

    \item \code{felm} now has some support for \code{"pdata.frame"}s from package plm.
    In particular, it is now possible to use things like \code{felm(y ~ x+
    lag(wage), data=pdata)} where \code{pdata} is a \code{"pdata.frame"}.

    \item A bug preventing the use of k-class estimators has been fixed.

    \item A bug leading to wrong robust standard errors when using
    weights has been fixed.
    
  }
}
\section{Changes in version 2.2 (2015-04-13)}{
  \itemize{

    \item A bug which caused \code{\link{felm}} to crash R when all
    cases of the data had NAs has been fixed. An error is raised instead.

    \item \code{\link{felm}} now supports weights.

    \item The conjugate gradient solver and trace estimation routine
    used by \code{\link{fevcov}} have been made public as
    \code{\link{cgsolve}} and \code{\link{mctrace}}.

    
  }
}

\section{Changes in version 2.0-2.1 (2015-03-05)}{
\itemize{
    \item Better support for matrices as covariates, both on the left
and right hand side of the equation, and in the instrumental variable
specification.

   \item Support for multiple left hand sides, with formulas like
\code{y|z|w ~ X |f1 + f2}.  The content of 'felm' objects have changed
slightly to accomodate this, and generic functions like vcov() and
coef() have gained an extra argument 'lhs' to specify which left hand
side to retrieve information about.

   \item Two new functions, waldtest() and condfstat() have been included to
test linear restrictions on parameters, and to compute conditional F
statistics to test for weak instruments in the IV estimations.
Observations with missing values are now correctly removed in IV-estimations.

   \item The old syntax, with 'G()' to specify factors to project out, has been
deprecated. It still works, but will yield a warning.  Similarly with
the arguments 'iv' and 'clustervar' to felm().

   \item \code{felm} now takes a new argument in the \code{...} part,
\code{kclass} to estimate k-class estimators with instrumental variables.
}
}

\section{Changes in version 1.8 (2014-12-25)}{
\itemize{
  \item Two new functions were added, fevcov() and bccorr(), to compute
limited mobility bias corrections for covariances (and correlations)
between high dimensional factors.  See citation('lfe') for reference
to article describing it.
}
}

\section{Changes in version 1.7 (2014-03-13)}{
\itemize{
   \item  Support for multiway clustering added.  Multiway clusters are specified
in the 4th part of the formula, like \code{y ~ x | f1+f2 | 0 | c1+c2}
}
}

\section{Changes in version 1.6 (2013-12-06)}{
\itemize{
   \item Added support for projecting out interactions between continuous
covariates and factors with specifications like \code{y ~ x | f1 + x2:f1}.

}
}

\section{Version 1.5 (2013-10-28)}{
\itemize{
   \item Version presented in R-journal article. (See citation('lfe') for
reference)
}
}
