# makematrix is a bit complicated. The purpose is to make model matrices for the various
# parts of the formulas.  The complications are due to the iv stuff.
# If there's an IV-part, its right hand side should be with the
# x. Their names are put in 'instruments'. Its left hand side goes in a separate entry 'ivy'



makematrix <- function(mf, contrasts = NULL, pf = parent.frame(),
                       clustervar = NULL, wildcard = "n", onlymm = FALSE) {
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  wpos <- which(!is.na(pmatch(names(mf), "weights")))
  if (length(wpos) > 0) {
    weights <- eval(mf[[wpos]], pf)
    if (!is.null(weights)) {
      if (anyNA(weights) || any(weights < 0)) stop("missing or negative weights not allowed")
      weights <- sqrt(weights)
      weights[weights == 0] <- 1e-60
    }
  } else {
    weights <- NULL
  }
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  # we should handle multiple lhs
  # but how?  model.frame() doesn't handle it, but we need
  # model.frame for subsetting and na.action, with the left hand side
  # included.  We create an artifical single lhs by summing the left hand
  # sides, just to get hold of the rhs.  Then we extract the left hand side

  # We need to remove the iv-spec from the Formula. It requires its own specification
  Form <- eval(mf[["formula"]], pf)
  formenv <- environment(Form)
  Form <- as.Formula(Form)

  # If Form rhs is shorter than 4, extend it with zeros.
  # Then we avoid some special cases later
  numrhs <- length(Form)[2]

  # we can't just dot-update the iv-part, update will only keep the instruments
  if (numrhs < 2) {
    Form <- update(Form, . ~ . | 0 | 0 | 0 | 0, drop = FALSE)
  } else if (numrhs < 3) {
    Form <- update(Form, . ~ . | . | 0 | 0 | 0, drop = FALSE)
  } else if (numrhs < 4) {
    # build from parts
    Form <- as.Formula(do.call(substitute, list(
      L ~ R1 | R2 | R3 | 0 | 0,
      list(
        L = formula(Form, lhs = NULL, rhs = 0)[[2]],
        R1 = formula(Form, lhs = 0, rhs = 1)[[2]],
        R2 = formula(Form, lhs = 0, rhs = 2)[[2]],
        R3 = formula(Form, lhs = 0, rhs = 3)[[2]]
      )
    )))
  } else if (numrhs < 5) {
    Form <- as.Formula(do.call(substitute, list(
      L ~ R1 | R2 | R3 | R4 | 0,
      list(
        L = formula(Form, lhs = NULL, rhs = 0)[[2]],
        R1 = formula(Form, lhs = 0, rhs = 1)[[2]],
        R2 = formula(Form, lhs = 0, rhs = 2)[[2]],
        R3 = formula(Form, lhs = 0, rhs = 3)[[2]],
        R4 = formula(Form, lhs = 0, rhs = 4)[[2]]
      )
    )))
  }


  if (numrhs > 5) stop("Formula can't have more than 5 parts")
  # Make a suitable formula for a model frame. No tricky IV-spec
  #  fullF <- formula(Form,lhs=NULL,rhs=0, drop=FALSE,collapse=TRUE,update=TRUE)
  fullF <- formula(Form, lhs = NULL, rhs = 0, drop = FALSE)
  for (i in seq_len(length(Form)[2])) {
    f <- formula(Form, lhs = 0, rhs = i, drop = FALSE)[[2]]
    if (i == 3) {
      if (identical(f, 0)) next
      f <- as.Formula(f[[2]]) # skip '('
      f <- formula(f, collapse = TRUE, drop = FALSE)
      fullF <- update(fullF, formula(substitute(. ~ . + F1 + F2, list(F1 = f[[2]], F2 = f[[3]]))), drop = FALSE, collapse = TRUE)
    } else {
      fullF <- update(fullF, formula(substitute(. ~ . + F, list(F = f))), drop = FALSE)
    }
  }

  usewild <- !identical(wildcard, "n")
  dataenv <- new.env(parent = pf)
  if (usewild) {
    # we must evalaute the data argument, but we want to
    # avoid it being reevaluated when we eval(mf),
    # so put it in an environment.  We do it like this
    # to have a short name in mf[['data']] in case of errors.
    data <- eval(mf[["data"]], pf)
    assign("..(@DATA@)..", data, dataenv)
    mf[["data"]] <- as.name("..(@DATA@)..")
    wildnames <- colnames(data)
    rm(data)
    if (wildcard == "R" || wildcard == "G") {
      wildnames <- unique(c(wildnames, rls(formenv)))
    }
    rewild <- wildcard %in% c("r", "R")
    fullF <- wildcard(fullF, wildnames, re = rewild)
  }
  environment(fullF) <- formenv
  mf[["formula"]] <- fullF

  # coerce pdata.frame (from plm) to ensure classes and attributes are preserved in model.frame
  # http://stackoverflow.com/questions/29724813/how-to-calculate-dynamic-panel-models-with-lfe-package
  if (!is.null(mf[["data"]])) {
    frname <- deparse(mf[["data"]])
    assign(
      "..pdata.coerce..",
      function(x) {
        if (inherits(x, "pdata.frame")) {
          if (!requireNamespace("plm")) {
            stop("Needs package plm to handle pdata.frame ", frname, call. = FALSE)
          }
          as.data.frame(x)
        } else {
          x
        }
      },
      dataenv
    )
    mf[["data"]] <- bquote(..pdata.coerce..(.(mf[["data"]])))
  }

  mfcall <- bquote(evalq(.(mf), .(dataenv)))
  mf <- eval(mfcall)

  if (nrow(mf) == 0) stop("0 (non-NA) cases; no valid data")

  rm(dataenv)
  naact <- na.action(mf)
  if (!is.null(naact) && !is.null(weights)) weights <- weights[-naact]

  #  if(is.null(mf$data)) data <- environment(mf[['formula']])
  # the factor list (rhs=2) needs special attention
  # it should be made into a model matrix, but treated specially.
  # It's a sum of terms like f + x:g

  fpart <- formula(Form, lhs = 0, rhs = 2)
  if (usewild) fpart <- wildcard(fpart, wildnames, re = rewild)
  ftm <- terms(fpart)

  # we make it into a call like
  # list(f=f, `x:g` = structure(g,x=x))
  # which we evaluate in the frame
  # make a function for ':'
  env <- new.env(parent = formenv)
  # make  '*' a function of two arguments to do the interaction.
  #  assign(':', function(a,b) {
  #    anam <- deparse(substitute(a))
  #    bnam <- deparse(substitute(b))
  #    message(' call : ',anam, ':' ,bnam)
  #    if(is.factor(a) && is.factor(b)) ret <- structure(interaction(a,b,drop=TRUE),xnam=bnam,fnam=anam)
  #    else if(is.factor(b)) ret <- structure(factor(b),x=a,xnam=anam,fnam=bnam)
  #    else if(is.factor(a)) ret <- structure(factor(a),x=b,xnam=bnam,fnam=anam)
  #    else stop('Error in term ',anam,':',bnam,'. Neither ',anam, ' nor ',bnam,' is a factor')
  #    ret
  #  }, env)
  #  fl <- eval(attr(ftm,'variables'), mf, env)
  vmat <- attr(ftm, "factors")

  fl <- lapply(attr(ftm, "term.labels"), function(tm) {
    # function for finding a factor name in the model frame.
    # It's really just to do mf[[n]], but in case of non-syntactical names like `a+b`,
    # the index name in mf is "a+b", whereas it's "`a+b'" in the terms object
    # so we must remove backticks before trying.
    gv <- function(n) mf[[sub("^`(.*)`$", "\\1", n)]]
    f <- gv(tm)

    # if it's a variable and only occurs in one term, pass it on
    if (!is.null(f) && sum(vmat[tm, ] > 0) == 1) {
      return(structure(factor(f), fnam = tm))
    }
    #    if(!is.null(f)) return(structure(factor(f),fnam=tm))
    # It's an interaction of some sort, find the variables in the interaction
    vars <- attr(ftm, "factors")[, tm]
    vars <- vars[vars != 0]
    nm <- names(vars)
    # find the factors
    isfac <- sapply(nm, function(n) is.factor(gv(n)))
    xx <- names(vars)[which(!isfac)]
    if (length(xx) > 1) stop("Interaction only allowed for one non-factor")
    # interact the factors
    # remove a reference level from the ones which are 1
    hasref <- vars == 1
    noref <- vars == 2

    # find the reference, we choose the largest level
    reffac <- which(isfac & hasref)
    namref <- names(vars[reffac])
    reflev <- sapply(namref, function(n) which(names(which.max(table(gv(n)))) %in% levels(gv(n))))
    names(reflev) <- namref
    # make a list with the reference level replaced
    if (length(xx) == 0) {
      #      rflist <- lapply(namref, function(n) {f <- mf[[n]]; levels(f)[[1]] <- NA; f})
      if (length(namref) == 1 && sum(noref & isfac) == 0) {
        rflist <- list(gv(namref))
      } else {
        ## GRM: Precursor to larger change below (i.e. don't replace individual
        ## FEs with NA yet, lest it creates too many reference cases once we
        ## interact them).
        # rflist <- lapply(namref, function(n) {f <- gv(n); levels(f)[[reflev[[n]]]] <- NA; f})
        rflist <- lapply(namref, function(n) {
          f <- gv(n)
          f
        })
      }
      names(rflist) <- namref
      ## GRM: Changed this next section of code to account for single reference
      ## case in the case of interacted FEs
      if (length(rflist) == 1) {
        f <- addNA(
          do.call(
            interaction,
            c(rflist, lapply(
              names(vars[noref & isfac]),
              function(n) gv(n)
            ), drop = TRUE)
          ),
          ifany = TRUE
        )
      } else {
        f <- do.call(interaction, c(rflist, lapply(
          names(vars[noref & isfac]),
          function(n) gv(n)
        ),
        drop = TRUE
        ))
        reflevcomb <- paste(reflev, collapse = ".")
        levels(f)[which(f == reflevcomb)] <- NA
        f <- addNA(f, ifany = TRUE)
      }
      refnam <- paste(sapply(namref, function(n) levels(gv(n))[reflev[n]]), collapse = "+")
      levels(f)[is.na(levels(f))] <- refnam
      #      structure(f, fnam=names(vars)[1], xnam=paste(names(vars)[-1],collapse=':'))
      structure(f, fnam = paste(names(vars), collapse = ":"))
    } else {
      f <- do.call(interaction, c(mf[names(vars)[isfac]], drop = TRUE))
      f <- f[!is.na(f)]
      structure(f, fnam = paste(names(vars[isfac]), collapse = ":"), x = mf[[xx]], xnam = xx)
    }
    #    f <- eval(parse(text=tm), mf, env)
    #    if(is.null(attr(f, 'fnam'))) factor(f) else f
  })
  names(fl) <- attr(ftm, "term.labels")

  # Name the interactions with the matrix first, then the factor name
  names(fl) <- sapply(names(fl), function(n) {
    f <- fl[[n]]
    x <- attr(f, "x", exact = TRUE)
    if (is.null(x)) {
      return(n)
    }
    return(paste(attr(f, "xnam"), attr(f, "fnam"), sep = ":"))
  })

  hasicpt <- all(sapply(fl, function(f) !is.null(attr(f, "x"))))
  environment(Form) <- formenv
  if (is.null(clustervar)) {
    cluform <- terms(formula(Form, lhs = 0, rhs = 4))
    cluster <- lapply(eval(attr(cluform, "variables"), mf, pf), factor)
    names(cluster) <- attr(cluform, "term.labels")
    if (length(cluster) == 0) cluster <- NULL
  } else {
    # backwards compatible
    if (is.character(clustervar)) clustervar <- as.list(clustervar)
    if (!is.list(clustervar)) clustervar <- list(clustervar)
    cluster <- lapply(clustervar, function(cv) {
      if (!is.character(cv)) factor(cv) else factor(eval(as.name(cv), mf, formenv))
    })
  }

  ivform <- formula(Form, lhs = 0, rhs = 3, drop = FALSE)

  # Pick up IV instruments
  if (ivform[[1]] == as.name("~")) ivform <- ivform[[2]]
  if (ivform[[1]] == as.name("(")) ivform <- ivform[[2]]
  if (!identical(ivform, 0)) {
    ivform <- as.Formula(ivform)
    if (length(ivform)[2] > 1) stop("Right hand side of IV-spec can't have multiple parts")
    inames <- as.character(attr(terms(formula(ivform, lhs = 0, rhs = 1)), "variables"))[-1]
    environment(ivform) <- formenv
  } else {
    ivform <- NULL
    inames <- NULL
  }

  # then the fifth part, the controls
  form <- formula(Form, lhs = 0, rhs = 5, drop = TRUE)
  if (!identical(form[[2]], 0)) {
    # always parse with intercept, remove it from matrix, so we never project out the intercept
    form <- formula(update(form, ~ . + 1))
    if (usewild) form <- wildcard(form, wildnames, re = rewild)
    ctrlterms <- terms(form, data = mf)
    ctrl <- delete.icpt(model.matrix(ctrlterms, data = mf, contrasts.arg = contrasts))
    if (typeof(ctrl) != "double") storage.mode(ctrl) <- "double"
    if (ncol(ctrl) == 0) {
      ctrlnames <- ctrl <- NULL
    } else {
      ctrlnames <- colnames(ctrl)
    }
  } else {
    ctrl <- NULL
    ctrlnames <- NULL
  }

  # We have taken Form apart. Keep only exogenous variables
  Form <- formula(Form, lhs = NULL, rhs = 1, drop = FALSE)
  environment(Form) <- formenv

  # model.response doesn't work with multiple responses
  #  y <- model.response(mf,"numeric")

  form <- formula(Form, lhs = NULL, rhs = 0, drop = FALSE)
  if (usewild) form <- wildcard(form, wildnames, re = rewild)
  y <- as.matrix(model.part(form, mf, lhs = NULL, rhs = 0), rownames.force = FALSE)
  if (typeof(y) != "double") storage.mode(y) <- "double"
  form <- formula(Form, lhs = 0, rhs = 1, collapse = c(FALSE, TRUE))
  if (usewild) form <- wildcard(form, wildnames, re = rewild)
  xterms <- terms(form, data = mf)
  x <- model.matrix(xterms, data = mf, contrasts.arg = contrasts)
  #  if(length(fl) > 0) {
  if (!hasicpt) {
    x <- delete.icpt(x)
    icpt <- FALSE
  } else {
    icpt <- attr(xterms, "intercept") != 0
  }
  if (typeof(x) != "double") storage.mode(x) <- "double"
  setdimnames(x, list(NULL, colnames(x)))

  if (!is.null(ivform)) {
    form <- formula(ivform, lhs = NULL, rhs = 0, drop = FALSE)
    if (usewild) form <- wildcard(form, wildnames, re = rewild)
    ivy <- as.matrix(model.part(form, mf, lhs = NULL, rhs = 0), rownames.force = FALSE)
    if (typeof(ivy) != "double") storage.mode(ivy) <- "double"
    form <- formula(ivform, lhs = 0, rhs = 1, collapse = c(FALSE, TRUE))
    if (usewild) form <- wildcard(form, wildnames, re = rewild)
    ivxterms <- terms(form, data = mf)
    # ivx should never contain an intercept
    ivx <- delete.icpt(model.matrix(ivxterms, data = mf, contrasts.arg = contrasts))
    if (typeof(ivx) != "double") storage.mode(ivx) <- "double"
    setdimnames(ivx, list(NULL, colnames(ivx)))
  } else {
    ivy <- NULL
    ivx <- NULL
  }

  mm <- list(x = x, y = y, ivx = ivx, ivy = ivy, ctrl = ctrl, fl = fl, weights = weights)
  mm$extra <- list(
    icpt = icpt, xterms = xterms, cluster = cluster, Form = Form, ivform = ivform,
    inames = inames, naact = naact, model = mf, mfcall = mfcall
  )
  if (onlymm) {
    return(mm)
  }
  mmdemean(mm)
}

mmdemean <- function(mm) {
  # orig is necessary to compute the r.residuals, i.e. residuals without dummies
  # it's used in getfe() and btrap, but is of no use if we have ctrl variables
  if (is.null(mm$weights)) {
    TSS <- apply(mm$y, 2, var) * (nrow(mm$y) - 1)
  } else {
    TSS <- apply(mm$y, 2, function(yy) sum(mm$weights^2 * (yy - sum(mm$weights^2 * yy / sum(mm$weights^2)))^2))
  }
  names(TSS) <- colnames(mm$y)

  if (length(mm$fl) != 0) {
    result <- demeanlist(list(y = mm$y, x = mm$x, ivy = mm$ivy, ivx = mm$ivx, ctrl = mm$ctrl),
      fl = mm$fl, weights = mm$weights
    )
    if (is.null(mm$ctrl)) result$orig <- list(y = mm$y, x = mm$x, ivy = mm$ivy, ivx = mm$ivx)
  } else {
    result <- list(y = mm$y, x = mm$x, ivy = mm$ivy, ivx = mm$ivx, ctrl = mm$ctrl)
  }

  if (!is.null(result$ctrl)) {
    # pure control variables to project out
    # do ols, use the residuals as new variables
    y <- cbind(result$y, result$x, result$ivy, result$ivx)
    x <- result$ctrl
    result$ctrl <- NULL
    #    fit <- .lm.fit(x,y)
    # my own is much faster for large datasets
    fit <- newols(list(y = y, x = x, weights = mm$weights), nostats = TRUE)
    resid <- as.matrix(fit$residuals)
    setdimnames(resid, list(NULL, colnames(y)))
    numctrl <- fit$rank
    rm(fit, x, y)

    result$y <- resid[, colnames(result$y), drop = FALSE]
    if (!is.null(result$x)) result$x <- resid[, colnames(result$x), drop = FALSE]
    if (!is.null(result$ivy)) result$ivy <- resid[, colnames(result$ivy), drop = FALSE]
    if (!is.null(result$ivx)) result$ivx <- resid[, colnames(result$ivx), drop = FALSE]
    rm(resid)
  } else {
    numctrl <- 0L
  }

  result$TSS <- TSS
  result$hasicpt <- mm$extra$icpt
  result$numctrl <- numctrl
  result$ctrlnames <- colnames(mm$ctrl)
  result$fl <- mm$fl
  result$terms <- mm$extra$xterms
  result$cluster <- mm$extra$cluster
  result$formula <- mm$extra$Form
  result$ivform <- mm$extra$ivform
  result$inames <- mm$extra$inames
  result$na.action <- mm$extra$naact
  result$weights <- mm$weights
  result$model <- mm$extra$model
  result$mfcall <- mm$extra$mfcall
  result
}


## Simple function borrowing from lme4::isNested() to check for nested factors.
## Will be used to check if a DoF correction needs to be made in the case where
## clusters are nested in FEs. See https://www.kellogg.northwestern.edu/faculty/matsa/htm/fe.htm
is_nested <- function(f1, f2) {
  f1 <- as.factor(f1)
  f2 <- as.factor(f2)
  stopifnot(length(f1) == length(f2))
  k <- length(levels(f1))
  sm <- as(methods::new("ngTMatrix", i = as.integer(f2) - 1L, j = as.integer(f1) -
    1L, Dim = c(length(levels(f2)), k)), "CsparseMatrix")
  all(sm@p[2:(k + 1L)] - sm@p[1:k] <= 1L)
}


newols <- function(mm, stage1 = NULL, pf = parent.frame(), nostats = FALSE, exactDOF = FALSE,
                   kappa = NULL, onlyse = FALSE, psdef = FALSE) {
  if (!is.null(mm$orig)) {
    orig <- mm$orig
  } else {
    orig <- mm
  }

  weights <- mm$weights

  numctrl <- if (is.null(mm$numctrl)) 0 else mm$numctrl
  hasicpt <- if (is.null(mm$hasicpt)) FALSE else mm$hasicpt
  cfactor <- compfactor(mm$fl)
  if (is.numeric(exactDOF)) {
    df <- exactDOF
    totvar <- nrow(mm$y) - df
  } else {
    # numrefs is also used later
    numrefs <- nrefs(mm$fl, cfactor, exactDOF)
    totvar <- totalpvar(mm$fl) - numrefs + numctrl
    df <- nrow(mm$y) - totvar
  }

  # special case for no covariates

  if (is.null(mm$x) || ncol(mm$x) == 0) {
    z <- list(
      N = nrow(mm$x),
      p = totvar, Pp = 0,
      na.action = mm$na.action, contrasts = mm$contrasts,
      df = df,
      nostats = FALSE,
      numctrl = numctrl,
      hasicpt = hasicpt,
      lhs = colnames(mm$y),
      call = match.call()
    )
    if (!onlyse) {
      z$r.residuals <- orig$y
      z$fe <- mm$fl
      z$cfactor <- cfactor
      z$fitted.values <- orig$y[, colnames(mm$y), drop = FALSE] - mm$y
      z$df.residual <- z$df
      z$residuals <- mm$y
      z$clustervar <- mm$cluster
    }
    class(z) <- "felm"
    return(z)
  }

  # lm.fit is an alternative.  Let's have a look at it later (didn't work before, but things have changed)
  # roll our own

  # to implement a k-class estimator, we should not project with P_Z, i.e.
  # onto the instruments. I.e. not X' (I-M_Z) X, but X' (I - kappa M_Z) X.
  # Indeed, the estimator is (X' (I-kappa M_Z)X)^{-1} X' (I-kappa M_Z) y)
  # Now, note that I - kappa M_Z = P_Z + (1-kappa)M_Z. So it is the
  # fitted values plus a fraction of the residuals

  # (see http://www.tandfonline.com/doi/pdf/10.1080/07350015.2014.978175 p 11)


  if (!is.null(weights)) iweights <- 1 / weights
  if (!is.null(weights)) {
    .Call(C_scalecols, mm$x, weights)
    .Call(C_scalecols, mm$y, weights)
  }

  if (!is.null(kappa)) {
    cp <- crossprod(mm$x) - kappa * crossprod(mm$noinst)
    b <- crossprod(mm$x, mm$y) - kappa * crossprod(mm$noinst, mm$y)
  } else {
    cp <- crossprod(mm$x)
    b <- crossprod(mm$x, mm$y)
  }
  ch <- cholx(cp)
  badvars <- attr(ch, "badvars")
  z <- list()
  class(z) <- "felm"
  if (is.null(badvars)) {
    beta <- backsolve(ch, backsolve(ch, b, transpose = TRUE))
    if (!nostats) z$inv <- chol2inv(ch)
  } else {
    beta <- matrix(NaN, nrow(cp), ncol(b))
    beta[-badvars, ] <- backsolve(ch, backsolve(ch, b[-badvars, ], transpose = TRUE))
    if (!nostats) {
      z$inv <- matrix(NA, nrow(cp), ncol(cp))
      z$inv[-badvars, -badvars] <- chol2inv(ch)
    }
  }
  if (!nostats && !is.null(kappa)) {
    # In k-class with k!=0 and k!=1, the covariance matrix isn't simply the
    # inverse of cp.  This is so because
    # hatbeta - beta = (X' K X)^{1} X' K' epsilon
    # Even when epsilon is iid, we obtain
    # var(hatbeta-beta) = sigma^2 (X' K X)^{-1} X' K' K X (X' K X)^{-1}
    # and since K isn't a projection, we do not have K'K = K, so
    # we can't cancel out one of the (X' K X)^{-1}
    #    kinv <- z$inv %*% crossprod(mm$x - kappa*mm$noinst) %*% z$inv
    kinv <- .Call(C_sandwich, 1.0, z$inv, crossprod(mm$x - kappa * mm$noinst))
  }
  rm(ch, b, cp)
  #  rownames(beta) <- colnames(orig$x)
  rownames(beta) <- colnames(mm$x)

  if (!is.null(weights)) {
    .Call(C_scalecols, mm$x, iweights)
    .Call(C_scalecols, mm$y, iweights)
  }

  #  z$lhs <- colnames(beta) <- colnames(orig$y)
  z$lhs <- colnames(beta) <- colnames(mm$y)
  z$hasicpt <- hasicpt
  z$TSS <- mm$TSS
  z$kappa <- kappa

  if (is.null(weights)) {
    z$P.TSS <- apply(mm$y, 2, var) * (nrow(mm$y) - 1)
  } else {
    z$P.TSS <- apply(mm$y, 2, function(yy) sum(weights^2 * (yy - sum(weights^2 * yy / sum(weights^2)))^2))
  }

  names(z$P.TSS) <- colnames(mm$y)
  if (!onlyse) z$weights <- weights

  z$numctrl <- numctrl
  z$coefficients <- z$beta <- beta
  # what else is there to put into a felm object?
  z$Pp <- ncol(orig$x)
  z$N <- nrow(orig$x)
  z$p <- z$Pp - length(badvars) + numctrl
  nabeta <- nazero(beta)

  zfit <- mm$x %*% nabeta
  zresid <- mm$y - zfit

  z$residuals <- zresid

  if (!onlyse) {
    z$response <- orig$y[, colnames(mm$y), drop = FALSE]
    z$c.fitted.values <- zfit
    z$fitted.values <- z$response - z$residuals
    #    z$fitted.values <- zfit
    z$cfactor <- compfactor(mm$fl)
    z$fe <- mm$fl
  }

  z$contrasts <- mm$contrasts
  if (!onlyse) {
    if (length(mm$fl) != 0) {
      #    message('dims:');print(dim(orig$y)); print(dim(orig$x)); print(dim(nabeta))
      if (is.null(kappa)) z$r.residuals <- orig$y - orig$x %*% nabeta
      #    if(!is.null(weights)) .Call(C_scalecols,z$r.residuals,weights)
    } else {
      z$r.residuals <- z$residuals
    }
  }

  # For IV, the residuals should be the residuals from the original
  # endogenous variables, not the predicted ones the difference are
  # the residuals from stage 1, which we must multiply by beta and
  # subtract.  the residuals from the 2nd stage are in iv.residuals
  # hmm, what about the r.residuals?  We modify them as well. They are
  # used in kaczmarz().

  if (!is.null(stage1)) {
    # we need the centred response in condfstat()
    fitnam <- makefitnames(stage1$lhs)
    ivresid <- stage1$residuals %*% nabeta[fitnam, , drop = FALSE]

    z$residuals <- z$residuals - ivresid
    if (!onlyse) {
      z$c.response <- mm$y
      z$iv.residuals <- zresid
      z$r.iv.residuals <- z$r.residuals
      z$r.residuals <- z$r.residuals - ivresid
      z$endovars <- stage1$lhs
      z$fitted.values <- z$response - z$residuals
    }
  }

  z$terms <- mm$terms
  totlev <- totalpvar(mm$fl)

  if (is.numeric(exactDOF)) {
    z$df <- exactDOF
    numdum <- z$N - z$p - z$df
    z$numrefs <- totlev - numdum
  } else {
    numdum <- totlev - numrefs
    z$numrefs <- numrefs
    z$df <- z$N - z$p - numdum
  }
  z$df.residual <- z$df
  z$rank <- z$N - z$df
  z$exactDOF <- exactDOF

  # should we subtract 1 for an intercept?
  # a similar adjustment is done in summary.felm when computing rdf
  z$p <- z$p + numdum #- 1
  z$xp <- z$p
  z$na.action <- mm$na.action
  class(z) <- "felm"
  cluster <- mm$cluster
  if (!onlyse) z$clustervar <- cluster
  z$stage1 <- stage1
  if (nostats) {
    z$nostats <- TRUE
    return(z)
  }

  z$nostats <- FALSE

  # then we go about creating the covariance matrices and tests
  # if there is a single lhs, they are just stored as matrices etc
  # in z.  If there are multiple lhs, these quantities are inserted
  # in a list z$STATS indexed by z$lhs
  # indexed by the name of the lhs

  vcvnames <- list(rownames(beta), rownames(beta))
  Ncoef <- nrow(beta)

  singlelhs <- length(z$lhs) == 1
  # preallocate STATS
  if (!singlelhs) z$STATS <- list()
  z$STATS <- list()
  if (is.null(kappa)) {
    vinv <- z$inv
  } else {
    vinv <- kinv
  }
  inv <- nazero(vinv)

  xz <- mm$x
  if (!is.null(kappa)) xz <- xz - kappa * mm$noinst
  for (lhs in z$lhs) {
    res <- z$residuals[, lhs]

    if (!is.null(weights)) res <- weights * res
    # when multiple lhs, vcvfactor is a vector
    # we need a list of vcvs in this case

    vcv <- sum(res**2) / z$df * vinv
    setdimnames(vcv, vcvnames)
    z$STATS[[lhs]]$vcv <- vcv
    if (singlelhs) z$vcv <- vcv

    # We should make the robust covariance matrix too.
    # it's inv * sum (X_i' u_i u_i' X_i) * inv
    # where u_i are the (full) residuals (Wooldridge, 10.5.4 (10.59))
    # i.e. inv * sum(u_i^2 X_i' X_i) * inv
    # for large datasets the sum is probably best computed by a series of scaled
    # rank k updates, i.e. the dsyrk blas routine, we make an R-version of it.
    # need to check this computation, the SE's are slightly numerically different from Stata's.
    # it seems stata does not do the small-sample adjustment
    dfadj <- z$N / z$df

    # Now, here's an optimzation for very large xz. If we use the wcrossprod and ccrossprod
    # functions, we can't get rid of xz, we end up with a copy of it which blows away memory.
    # we need to scale xz with the residuals in xz, but we don't want to expand res to a full matrix,
    # and even get a copy in the result.
    # Thus we modify it in place with a .Call. The scaled variant is also used in the cluster computation.

    rscale <- ifelse(res == 0, 1e-40, res) # make sure nothing is zero
    if (!is.null(weights)) rscale <- rscale * weights
    # This one scales the columns without copying
    # For xz, remember to scale it back, because we scale directly into
    # mm$x
    .Call(C_scalecols, xz, rscale)
    # compute inv %*% crossprod(xz) %*% inv
    # via a blas dsyrk. Save some memory
    meat <- matrix(0, Ncoef, Ncoef)
    .Call(C_dsyrk, 0.0, meat, dfadj, xz)
    rvcv <- .Call(C_sandwich, 1.0, inv, meat)
    setdimnames(rvcv, vcvnames)

    z$STATS[[lhs]]$robustvcv <- rvcv
    if (singlelhs) z$robustvcv <- rvcv

    rm(meat, rvcv)

    # then the clustered covariance matrix
    if (!is.null(cluster)) {
      method <- attr(cluster, "method")
      if (is.null(method)) method <- "cgm"
      dfadj <- (z$N - 1) / z$df

      ## GRM: Extra adjustments to the DoF are needed in cases where clusters
      ## are nested within any of the FEs. See Cameron and Miller (2015, pp. 14-15):
      ## http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf#page=14
      fe_cl_grid <- expand.grid(fe_k = seq_along(z$fe), cl_g = seq_along(cluster))
      any_nested <-
        vapply(
          seq_len(nrow(fe_cl_grid)),
          function(n) {
            fe_k <- fe_cl_grid$fe_k[n]
            cl_g <- fe_cl_grid$cl_g[n]
            is_nested(z$fe[[fe_k]], cluster[[cl_g]])
          },
          FUN.VALUE = logical(1)
        )
      ## Find the minimum cluster dimension. Will be used below in the case of
      ## multiway clustering, but only if the FEs are nested within a cluster,
      ## or 'cgm2' (or 'reghdfe') is specified for the `cmethod` argument.
      min_clust <-
        min(vapply(
          seq_along(cluster),
          function(i) nlevels(cluster[[i]]),
          FUN.VALUE = integer(1)
        ))
      if (any(any_nested)) {
        # Will use the simple correction proposed by Gormley and Matsa (RFS, 2014).
        # See: https://www.kellogg.northwestern.edu/faculty/matsa/htm/fe.htm
        dfadj <- dfadj * z$df / (z$df + totvar - 1)
        # In addition to the above, the nested clusters case requires that
        # regressor p-values are calculated using no. of clusters - 1 degrees
        # of freedom; similar to "df2" in `waldtest()`. This is straight forward
        # when there is only a single cluster variable. In the case of multiway
        # clustering, however, we'll conservatively take "no. of clusters" to
        # mean the cluster variable with the smallest dimension. If nothing else,
        # this should ensure consistency with comparable implementations in
        # Stata (via reghdfe) and Julia (via FixedEffectModels.jl). See also:
        # https://github.com/matthieugomez/FixedEffectModels.jl/pull/50
        z$df <- min(min_clust - 1, z$df)
        z$df.residual <- z$df
      }
      ## End of nested cluster adjustment

      d <- length(cluster)
      if (method %in% c("cgm", "cgm2", "reghdfe")) {
        meat <- matrix(0, Ncoef, Ncoef)
        for (i in 1:(2^d - 1)) {
          # Find out which ones to interact
          iac <- as.logical(intToBits(i))[1:d]
          # odd number is positive, even is negative
          sgn <- 2 * (sum(iac) %% 2) - 1
          # interact the factors
          ia <- factor(do.call(paste, c(cluster[iac], sep = "\004")))
          # CGM2011 (sec 2.3) describe two possible small-sample adjustments
          # using the number of clusters in each cluster category. Note that
          # these two approaches should only diverge in the case of multiway
          # clustering.
          if (method == "cgm") {
            ## Option 1 (used by GCM2011 in their paper and also our default here)
            adj <- sgn * dfadj * nlevels(ia) / (nlevels(ia) - 1)
          } else { ## i.e. if method %in% c('cgm2','reghdfe')
            ## Option 2 (used by Stata's reghdfe among others, so we'll give it that alias for convenience)
            adj <- sgn * dfadj * min_clust / (min_clust - 1)
            ## Will also need to adjust DoF used to calculate p-vals and CIs
            z$df <- min(min_clust - 1, z$df)
            z$df.residual <- z$df
          }
          .Call(C_dsyrk, 1.0, meat, adj, crowsum(xz, ia))
        }

        cvcv <- .Call(C_sandwich, 1.0, inv, meat)
        if (psdef && d > 1) {
          ev <- eigen(cvcv)
          badev <- Im(ev$values) != 0 | Re(ev$values) < 0
          if (any(badev)) {
            warning("Negative eigenvalues set to zero in multiway clustered variance matrix. See felm(...,psdef=FALSE)")
            ev$values[badev] <- 0
            cvcv <- Re(ev$vectors %*% diag(ev$values) %*% t(ev$vectors))
          }
          #          if(any(Im(ev$values) != 0 | Re(ev$values) < 0)) {
          #            warning('Negative eigenvalues set to zero in multiway clustered variance matrix. See felm(...,psdef=FALSE)')
          #            cvcv <- ev$vectors %*% diag(pmax(ev$values,0)) %*% t(ev$vectors)
          #          }
          rm(ev)
        }
        setdimnames(cvcv, vcvnames)
        z$STATS[[lhs]]$clustervcv <- cvcv
        if (singlelhs) z$clustervcv <- cvcv
        rm(meat, cvcv)
      } else if (method == "gaure") {
        #       .Call(C_scalecols, xz, 1/rscale)
        meat <- matrix(0, nrow(z$vcv), ncol(z$vcv))
        # scale the columns according to group size
        sc <- apply(sapply(cluster, function(f) table(f)[f]), 1, mean)
        xc <- demeanlist(xz, cluster, means = TRUE)
        .Call(C_scalecols, xc, sqrt(sc))
        adj <- dfadj
        #        adj <- adj*prod(sapply(cluster, function(f) nlevels(f)/(nlevels(f)-1)))
        .Call(C_dsyrk, 1, meat, adj, xc)
        cvcv <- .Call(C_sandwich, 1.0, inv, meat)
        setdimnames(cvcv, vcvnames)
        z$STATS[[lhs]]$clustervcv <- cvcv
        if (singlelhs) z$clustervcv <- cvcv
        rm(meat, cvcv)
        #        .Call(C_scalecols, xz, rscale)
      } else {
        stop("unknown multi way cluster algorithm:", method)
      }


      z$STATS[[lhs]]$cse <- sqrt(diag(z$STATS[[lhs]]$clustervcv))
      z$STATS[[lhs]]$ctval <- z$coefficients[, lhs] / z$STATS[[lhs]]$cse
      z$STATS[[lhs]]$cpval <- 2 * pt(abs(z$STATS[[lhs]]$ctval), z$df, lower.tail = FALSE)

      if (singlelhs) {
        z$cse <- z$STATS[[lhs]]$cse
        z$ctval <- z$STATS[[lhs]]$ctval
        z$cpval <- z$STATS[[lhs]]$cpval
      }
    }

    z$STATS[[lhs]]$se <- sqrt(diag(z$STATS[[lhs]]$vcv))
    z$STATS[[lhs]]$tval <- z$coefficients[, lhs] / z$STATS[[lhs]]$se
    z$STATS[[lhs]]$pval <- 2 * pt(abs(z$STATS[[lhs]]$tval), z$df, lower.tail = FALSE)

    z$STATS[[lhs]]$rse <- sqrt(diag(z$STATS[[lhs]]$robustvcv))
    z$STATS[[lhs]]$rtval <- z$coefficients[, lhs] / z$STATS[[lhs]]$rse
    z$STATS[[lhs]]$rpval <- 2 * pt(abs(z$STATS[[lhs]]$rtval), z$df, lower.tail = FALSE)

    if (singlelhs) {
      z$se <- z$STATS[[lhs]]$se
      z$tval <- z$STATS[[lhs]]$tval
      z$pval <- z$STATS[[lhs]]$pval

      z$rse <- z$STATS[[lhs]]$rse
      z$rtval <- z$STATS[[lhs]]$rtval
      z$rpval <- z$STATS[[lhs]]$rpval
    }

    # reset this for next lhs
    .Call(C_scalecols, xz, 1 / rscale)
  }
  if (onlyse) z$residuals <- NULL
  z
}







#' Fit a linear model with multiple group fixed effects
#'
#' 'felm' is used to fit linear models with multiple group fixed effects,
#' similarly to lm.  It uses the Method of Alternating projections to sweep out
#' multiple group effects from the normal equations before estimating the
#' remaining coefficients with OLS.
#'
#' This function is intended for use with large datasets with multiple group
#' effects of large cardinality.  If dummy-encoding the group effects results
#' in a manageable number of coefficients, you are probably better off by using
#' [lm()].
#'
#' The formula specification is a response variable followed by a four part
#' formula. The first part consists of ordinary covariates, the second part
#' consists of factors to be projected out. The third part is an
#' IV-specification. The fourth part is a cluster specification for the
#' standard errors.  I.e. something like `y ~ x1 + x2 | f1 + f2 | (Q|W ~
#' x3+x4) | clu1 + clu2` where `y` is the response, `x1,x2` are
#' ordinary covariates, `f1,f2` are factors to be projected out, `Q`
#' and `W` are covariates which are instrumented by `x3` and
#' `x4`, and `clu1,clu2` are factors to be used for computing cluster
#' robust standard errors. Parts that are not used should be specified as
#' `0`, except if it's at the end of the formula, where they can be
#' omitted.  The parentheses are needed in the third part since `|` has
#' higher precedence than `~`. Multiple left hand sides like `y|w|x ~
#' x1 + x2 |f1+f2|...` are allowed.
#'
#' Interactions between a covariate `x` and a factor `f` can be
#' projected out with the syntax `x:f`.  The terms in the second and
#' fourth parts are not treated as ordinary formulas, in particular it is not
#' possible with things like `y ~ x1 | x*f`, rather one would specify
#' `y ~ x1 + x | x:f + f`. Note that `f:x` also works, since R's
#' parser does not keep the order.  This means that in interactions, the factor
#' *must* be a factor, whereas a non-interacted factor will be coerced to
#' a factor. I.e. in `y ~ x1 | x:f1 + f2`, the `f1` must be a factor,
#' whereas it will work as expected if `f2` is an integer vector.
#'
#' In older versions of \pkg{lfe} the syntax was `felm(y ~ x1 + x2 + G(f1)
#' + G(f2), iv=list(Q ~ x3+x4, W ~ x3+x4), clustervar=c('clu1','clu2'))`. This
#' syntax still works, but yields a warning. Users are *strongly*
#' encouraged to change to the new multipart formula syntax.  The old syntax
#' will be removed at a later time.
#'
#' The standard errors are adjusted for the reduced degrees of freedom coming
#' from the dummies which are implicitly present.  (An exception occurs in the
#' case of clustered standard errors and, specifically, where clusters are
#' nested within fixed effects; see
#' [here](https://github.com/sgaure/lfe/issues/1#issuecomment-528643802).)
#' In the case of two factors,
#' the exact number of implicit dummies is easy to compute.  If there are more
#' factors, the number of dummies is estimated by assuming there's one
#' reference-level for each factor, this may be a slight over-estimation,
#' leading to slightly too large standard errors. Setting `exactDOF='rM'`
#' computes the exact degrees of freedom with `rankMatrix()` in package
#' \pkg{Matrix}.
#'
#' For the iv-part of the formula, it is only necessary to include the
#' instruments on the right hand side.  The other explanatory covariates, from
#' the first and second part of `formula`, are added automatically in the
#' first stage regression.  See the examples.
#'
#' The `contrasts` argument is similar to the one in `lm()`, it is
#' used for factors in the first part of the formula. The factors in the second
#' part are analyzed as part of a possible subsequent `getfe()` call.
#'
#' The `cmethod` argument may affect the clustered covariance matrix (and
#' thus regressor standard errors), either directly or via adjustments to a
#' degrees of freedom scaling factor. In particular, Cameron, Gelbach and Miller
#' (CGM2011, sec. 2.3) describe two possible small cluster corrections that are
#' relevant in the case of multiway clustering. \itemize{
#' \item The first approach adjusts each component of the cluster-robust
#' variance estimator (CRVE) by its own \eqn{c_i} adjustment factor. For
#' example, the first component (with \eqn{G} clusters) is adjusted by
#' \eqn{c_1=\frac{G}{G-1}\frac{N-1}{N-K}}{c_1 = G/(G-1)*(N-1)/(N-K)},
#' the second component (with \eqn{H} clusters) is adjusted
#' by \eqn{c_2=\frac{H}{H-1}\frac{N-1}{N-K}}{c_2 = H/(H-1)*(N-1)/(N-K)}, etc.
#' \item The second approach applies the same adjustment to all CRVE components:
#' \eqn{c=\frac{J}{J-1}\frac{N-1}{N-K}}{c = J/(J-1)*(N-1)/(N-K)}, where
#' \eqn{J=\min(G,H)}{J=min(G,H)} in the case of two-way clustering, for example.
#' }
#' Any differences resulting from these two approaches are likely to be minor,
#' and they will obviously yield exactly the same results when there is only one
#' cluster dimension. Still, CGM2011 adopt the former approach in their own
#' paper and simulations. This is also the default method that `felm` uses
#' (i.e. `cmethod = 'cgm'`). However, the latter approach has since been
#' adopted by several other packages that allow for robust inference with
#' multiway clustering. This includes the popular Stata package
#' [reghdfe](http://scorreia.com/software/reghdfe/), as well as the
#' [FixedEffectModels.jl](https://github.com/matthieugomez/FixedEffectModels.jl)
#' implementation in Julia. To match results from these packages exactly, use
#' `cmethod = 'cgm2'` (or its alias, `cmethod = 'reghdfe'`). It is
#' possible that some residual differences may still remain; see discussion
#' [here](https://github.com/sgaure/lfe/issues/1#issuecomment-530561314).
#'
#' The old syntax with a single part formula with the `G()` syntax for the
#' factors to transform away is still supported, as well as the
#' `clustervar` and `iv` arguments, but users are encouraged to move
#' to the new multi part formulas as described here.  The `clustervar` and
#' `iv` arguments have been moved to the `...` argument list. They
#' will be removed in some future update.
#'
#' @param formula an object of class '"formula"' (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted. Similarly to
#' 'lm'.  See Details.
#' @param data a data frame containing the variables of the model.
#' @param exactDOF logical. If more than two factors, the degrees of freedom
#' used to scale the covariance matrix (and the standard errors) is normally
#' estimated. Setting `exactDOF=TRUE` causes `felm` to attempt to
#' compute it, but this may fail if there are too many levels in the factors.
#' `exactDOF='rM'` will use the exact method in
#' `Matrix::rankMatrix()`, but this is slower. If neither of these methods
#' works, it is possible to specify `exactDOF='mc'`, which utilizes a
#' Monte-Carlo method to estimate the expectation E(x' P x) = tr(P), the trace
#' of a certain projection, a method which may be more accurate than the
#' default guess.
#'
#' If the degrees of freedom for some reason are known, they can be specified
#' like `exactDOF=342772`.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain `NA`s.  The default is set by the `na.action` setting of
#' `options`, and is `na.fail` if that is unset.  The 'factory-fresh'
#' default is `na.omit`.  Another possible value is `NULL`, no
#' action. `na.exclude` is currently not supported.
#' @param contrasts an optional list. See the `contrasts.arg` of
#' `model.matrix.default`.
#' @param weights an optional vector of weights to be used in the fitting
#' process.  Should be 'NULL' or a numeric vector.  If non-NULL, weighted least
#' squares is used with weights `weights` (that is, minimizing
#' `sum(w*e^2)`); otherwise ordinary least squares is used.
#' @param ... other arguments.  \itemize{
#'
#' \item `cmethod` character. Which clustering method to use. Known
#' arguments are `'cgm'` (the default), `'cgm2'` (or `'reghdfe'`,
#' its alias). These alternate methods will generally
#' yield equivalent results, except in the case of multiway clustering with few
#' clusters along at least one dimension.
#'
#' \item `keepX` logical. To include a copy of the expanded data matrix in
#' the return value, as needed by [bccorr()] and [fevcov()]
#' for proper limited mobility bias correction.
#'
#' \item `keepCX` logical. Keep a copy of the centred expanded data matrix
#' in the return value. As list elements `cX` for the explanatory
#' variables, and `cY` for the outcome.
#'
#' \item `keepModel` logical. Keep a copy of the model frame.
#'
#' \item `nostats` logical. Don't include covariance matrices in the
#' output, just the estimated coefficients and various descriptive information.
#' For IV, `nostats` can be a logical vector of length 2, with the last
#' value being used for the 1st stages.  \item `psdef` logical. In case of
#' multiway clustering, the method of Cameron, Gelbach and Miller may yield a
#' non-definite variance matrix. Ordinarily this is forced to be semidefinite
#' by setting negative eigenvalues to zero. Setting `psdef=FALSE` will
#' switch off this adjustment.  Since the variance estimator is asymptotically
#' correct, this should only have an effect when the clustering factors have
#' very few levels.
#'
#' \item `kclass` character. For use with instrumental variables. Use a
#' k-class estimator rather than 2SLS/IV. Currently, the values `'nagar',
#' 'b2sls', 'mb2sls', 'liml'` are accepted, where the names are from
#' \cite{Kolesar et al (2014)}, as well as a numeric value for the 'k' in
#' k-class. With `kclass='liml'`, `felm` also accepts the argument
#' `fuller=<numeric>`, for using a Fuller adjustment of the
#' liml-estimator.
#'
#' \item `Nboot, bootexpr, bootcluster` Since `felm` has quite a bit
#' of overhead in the creation of the model matrix, if one wants confidence
#' intervals for some function of the estimated parameters, it is possible to
#' bootstrap internally in `felm`.  That is, the model matrix is resampled
#' `Nboot` times and estimated, and the `bootexpr` is evaluated
#' inside an `sapply`. The estimated coefficients and the left hand
#' side(s) are available by name. Any right hand side variable `x` is
#' available by the name `var.x`.  The `"felm"`-object for each
#' estimation is available as `est`. If a `bootcluster` is specified
#' as a factor, entire levels are resampled. `bootcluster` can also be a
#' function with no arguments, it should return a vector of integers, the rows
#' to use in the sample. It can also be the string 'model', in which case the
#' cluster is taken from the model. `bootexpr` should be an expression,
#' e.g. like `quote(x/x2 * abs(x3)/mean(y))`.  It could be wise to specify
#' `nostats=TRUE` when bootstrapping, unless the covariance matrices are
#' needed in the bootstrap. If you need the covariance matrices in the full
#' estimate, but not in the bootstrap, you can specify it in an attribute
#' `"boot"` as `nostats=structure(FALSE, boot=TRUE)`.
#'
#' \item `iv, clustervar` deprecated.  These arguments will be removed at
#' a later time, but are still supported in this field. Users are
#' *STRONGLY* encouraged to use multipart formulas instead.  In
#' particular, not all functionality is supported with the deprecated syntax;
#' iv-estimations actually run a lot faster if multipart formulas are used, due
#' to new algorithms which I didn't bother to shoehorn in place for the
#' deprecated syntax.
#'
#' }
#' @return `felm` returns an object of `class` `"felm"`.  It is
#' quite similar to an `"lm"` object, but not entirely compatible.
#'
#' The generic `summary`-method will yield a summary which may be
#' `print`'ed.  The object has some resemblance to an `'lm'` object,
#' and some postprocessing methods designed for `lm` may happen to work.
#' It may however be necessary to coerce the object to succeed with this.
#'
#' The `"felm"` object is a list containing the following fields:
#'
#' \item{coefficients}{a numerical vector. The estimated coefficients.}
#' \item{N}{an integer. The number of observations} \item{p}{an integer. The
#' total number of coefficients, including those projected out.}
#' \item{response}{a numerical vector. The response vector.}
#' \item{fitted.values}{a numerical vector. The fitted values.}
#'
#' \item{residuals}{a numerical vector. The residuals of the full system, with
#' dummies.  For IV-estimations, this is the residuals when the original
#' endogenous variables are used, not their predictions from the 1st stage.}
#'
#' \item{r.residuals}{a numerical vector. Reduced residuals, i.e. the residuals
#' resulting from predicting *without* the dummies.}
#'
#' \item{iv.residuals}{numerical vector. When using instrumental variables,
#' residuals from 2. stage, i.e. when predicting with the predicted endogenous
#' variables from the 1st stage.}
#'
#' \item{weights}{numeric. The square root of the argument `weights`.}
#'
#' \item{cfactor}{factor of length N. The factor describing the connected
#' components of the two first terms in the second part of the model formula.}
#'
#' \item{vcv}{a matrix. The variance-covariance matrix.}
#'
#' \item{fe}{list of factors. A list of the terms in the second part of the
#' model formula.}
#'
#' \item{stage1}{The '`felm`' objects for the IV 1st stage, if used. The
#' 1st stage has multiple left hand sides if there are more than one
#' instrumented variable.}
#'
#' \item{iv1fstat}{list of numerical vectors. For IV 1st stage, F-value for
#' excluded instruments, the number of parameters in restricted model and in
#' the unrestricted model.}
#'
#' \item{X}{matrix. The expanded data matrix, i.e. from the first part of the
#' formula. To save memory with large datasets, it is only included if
#' `felm(keepX=TRUE)` is specified.  Must be included if
#' [bccorr()] or [fevcov()] is to be used for correcting
#' limited mobility bias.  }
#'
#' \item{cX, cY}{matrix. The centred expanded data matrix. Only included if
#' `felm(keepCX=TRUE)`.  }
#'
#' \item{boot}{The result of a `replicate` applied to the `bootexpr`
#' (if used).}
#'
#' @note
#' Side effect: If `data` is an object of class `"pdata.frame"` (from
#' the \pkg{plm} package), the \pkg{plm} namespace is loaded if available, and
#' `data` is coerced to a `"data.frame"` with `as.data.frame`
#' which dispatches to a \pkg{plm} method.  This ensures that transformations
#' like `diff` and `lag` from \pkg{plm} works as expected, but it
#' also incurs an additional copy of the `data`, and the \pkg{plm}
#' namespace remains loaded after `felm` returns.  When working with
#' `"pdata.frame"`s, this is what is usually wanted anyway.
#'
#' For technical reasons, when running IV-estimations, the data frame supplied
#' in the `data` argument to `felm`, should *not* contain
#' variables with names ending in `'(fit)'`.  Variables with such names
#' are used internally by `felm`, and may then accidentally be looked up
#' in the data frame instead of the local environment where they are defined.
#' @seealso [getfe()] [summary.felm()]
#' [condfstat()] [waldtest()]
#' @references Cameron, A.C., J.B. Gelbach and D.L. Miller (2011) \cite{Robust
#' inference with multiway clustering}, Journal of Business & Economic
#' Statistics 29 (2011), no. 2, 238--249.
#' \doi{10.1198/jbes.2010.07136}
#'
#' Kolesar, M., R. Chetty, J. Friedman, E. Glaeser, and G.W. Imbens (2014)
#' \cite{Identification and Inference with Many Invalid Instruments}, Journal
#' of Business & Economic Statistics (to appear).
#' \doi{10.1080/07350015.2014.978175}
#' @examples
#'
#' ## Default is to use all cores. We'll limit it to 2 for this example.
#' oldopts <- options("lfe.threads")
#' options(lfe.threads = 2)
#'
#' ## Simulate data
#' set.seed(42)
#' n <- 1e3
#'
#' d <- data.frame(
#'   # Covariates
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   # Individuals and firms
#'   id = factor(sample(20, n, replace = TRUE)),
#'   firm = factor(sample(13, n, replace = TRUE)),
#'   # Noise
#'   u = rnorm(n)
#' )
#'
#' # Effects for individuals and firms
#' id.eff <- rnorm(nlevels(d$id))
#' firm.eff <- rnorm(nlevels(d$firm))
#'
#' # Left hand side
#' d$y <- d$x1 + 0.5 * d$x2 + id.eff[d$id] + firm.eff[d$firm] + d$u
#'
#' ## Estimate the model and print the results
#' est <- felm(y ~ x1 + x2 | id + firm, data = d)
#' summary(est)
#' # Compare with lm
#' summary(lm(y ~ x1 + x2 + id + firm - 1, data = d))
#'
#' ## Example with 'reverse causation' (IV regression)
#'
#' # Q and W are instrumented by x3 and the factor x4.
#' d$x3 <- rnorm(n)
#' d$x4 <- sample(12, n, replace = TRUE)
#' d$Q <- 0.3 * d$x3 + d$x1 + 0.2 * d$x2 + id.eff[d$id] + 0.3 * log(d$x4) - 0.3 * d$y +
#'   rnorm(n, sd = 0.3)
#' d$W <- 0.7 * d$x3 - 2 * d$x1 + 0.1 * d$x2 - 0.7 * id.eff[d$id] + 0.8 * cos(d$x4) -
#'   0.2 * d$y + rnorm(n, sd = 0.6)
#'
#' # Add them to the outcome variable
#' d$y <- d$y + d$Q + d$W
#'
#' ## Estimate the IV model and report robust SEs
#' ivest <- felm(y ~ x1 + x2 | id + firm | (Q | W ~ x3 + factor(x4)), data = d)
#' summary(ivest, robust = TRUE)
#' condfstat(ivest)
#' # Compare with the not instrumented fit:
#' summary(felm(y ~ x1 + x2 + Q + W | id + firm, data = d))
#'
#' ## Example with multiway clustering
#'
#' # Create a large cluster group (500 clusters) and a small one (20 clusters)
#' d$cl1 <- factor(sample(rep(1:500, length.out = n)))
#' d$cl2 <- factor(sample(rep(1:20, length.out = n)))
#' # Function for adding clustered noise to our outcome variable
#' cl_noise <- function(cl) {
#'   obs_per_cluster <- n / nlevels(cl)
#'   unlist(replicate(nlevels(cl),
#'     rnorm(obs_per_cluster, mean = rnorm(1), sd = runif(1)),
#'     simplify = FALSE
#'   ))
#' }
#'
#' # New outcome variable
#' d$y_cl <- d$x1 + 0.5 * d$x2 + id.eff[d$id] + firm.eff[d$firm] +
#'   cl_noise(d$cl1) + cl_noise(d$cl2)
#'
#' ## Estimate and print the model with cluster-robust SEs (default)
#' est_cl <- felm(y_cl ~ x1 + x2 | id + firm | 0 | cl1 + cl2, data = d)
#' summary(est_cl)
#'
#' # Print ordinary standard errors:
#' summary(est_cl, robust = FALSE)
#' # Match cluster-robust SEs from Stata's reghdfe package:
#' summary(felm(y_cl ~ x1 + x2 | id + firm | 0 | cl1 + cl2,
#'   data = d,
#'   cmethod = "reghdfe"
#' ))
#'
#' ## Restore default options
#' options(oldopts)
#'
#' @export felm
felm <- function(formula, data, exactDOF = FALSE, subset, na.action,
                 contrasts = NULL, weights = NULL, ...) {
  knownargs <- c(
    "iv", "clustervar", "cmethod", "keepX", "nostats",
    "wildcard", "kclass", "fuller", "keepCX", "Nboot",
    "bootexpr", "bootcluster", "onlyse", "psdef", "keepModel"
  )
  cmethod <- "cgm"
  iv <- NULL
  clustervar <- NULL
  nostats <- FALSE
  wildcard <- "n"
  kclass <- NULL
  fuller <- 0
  Nboot <- 0
  onlyse <- FALSE
  bootexpr <- NULL
  bootcluster <- NULL
  deprec <- c("iv", "clustervar")
  psdef <- TRUE
  keepX <- FALSE
  keepCX <- FALSE
  keepModel <- FALSE

  mf <- match.call(expand.dots = TRUE)

  # Currently there shouldn't be any ... arguments
  # check that the list is empty

  #  if(length(mf[['...']]) > 0) stop('unknown argument ',mf['...'])


  args <- list(...)
  ka <- knownargs[pmatch(names(args), knownargs, duplicates.ok = FALSE)]
  names(args)[!is.na(ka)] <- ka[!is.na(ka)]
  dpr <- deprec[match(ka, deprec)]
  if (any(!is.na(dpr))) {
    bad <- dpr[which(!is.na(dpr))]
    .Deprecated("", msg = paste("Argument(s)", paste(bad, collapse = ","), "are deprecated and will be removed, use multipart formula instead"))
    #    warning('Argument(s) ',paste(bad,collapse=','), ' are deprecated and will be removed, use multipart formula instead')
  }
  env <- environment()
  lapply(intersect(knownargs, ka), function(arg) assign(arg, args[[arg]], pos = env))

  if (!(cmethod %in% c("cgm", "cgm2", "reghdfe", "gaure"))) stop("Unknown cmethod: ", cmethod)

  # also implement a check for unknown arguments
  unk <- setdiff(names(args), knownargs)
  if (length(unk) > 0) stop("unknown arguments ", paste(unk, collapse = " "))

  # backwards compatible
  Gtm <- terms(formula(as.Formula(formula), rhs = 1), specials = "G")
  if (!is.null(attr(Gtm, "specials")$G) || !is.null(iv)) {
    mf <- match.call(expand.dots = TRUE)
    mf[[1L]] <- quote(..oldfelm)
    return(eval.parent(mf))
  }
  pint <- getOption("lfe.pint")
  start <- last <- Sys.time()
  mm <- makematrix(mf, contrasts, pf = parent.frame(), clustervar, wildcard = wildcard)
  if (!is.null(mm$cluster)) attr(mm$cluster, "method") <- cmethod
  now <- Sys.time()
  if (now > last + pint) {
    last <- now
    message(date(), " finished centering model matrix")
  }
  z <- felm.mm(mm, nostats, exactDOF, keepX, keepCX, keepModel, kclass, fuller, onlyse, psdef = psdef)
  z$call <- match.call()
  z$formula <- formula
  z$keepX <- keepX
  z$keepCX <- keepCX
  if (Nboot > 0) {
    now <- Sys.time()
    if (now > last + pint) {
      last <- now
      message(date(), " finished estimate, starting bootstrap")
    }

    mm <- makematrix(mf, contrasts,
      pf = parent.frame(), clustervar, wildcard = wildcard,
      onlymm = TRUE
    )
    if (is.null(bootexpr)) bootexpr <- quote(beta)
    if (is.null(bootcluster)) {
      csample <- function() sort(sample(nrow(mm$x), replace = TRUE))
    } else if (is.function(bootcluster)) {
      csample <- bootcluster
    } else if (is.factor(bootcluster)) {
      csample <- function() resample(bootcluster, na.action = mm$extra$naact)
    } else if (identical(bootcluster, "model")) {
      csample <- function() resample(mm$extra$cluster)
    } else {
      stop("bootcluster must be either a factor or a function")
    }
    bootstat <- nostats
    if (!is.null(attr(nostats, "boot"))) bootstat <- attr(nostats, "boot")
    iii <- 0
    bootfun <- function() {
      assign("now", Sys.time(), inherits = TRUE) # replaces: now <<- Sys.time()
      iii_tmp <- iii + 1 # replaces: iii <<- iii+1
      assign("iii", iii_tmp, inherits = TRUE)
      if (now > last + pint) {
        assign("last", now, inherits = TRUE) # replaces last <<- now
        message(date(), " done boot iter ", iii)
      }
      bootenv <- new.env()
      # we delay assign to avoid unnecessary estimating and copying
      if (FALSE) {
        olsmms <- mms <- est <- bootx <- booty <- bootivy <- NULL
      } # avoid warning about no visible binding
      delayedAssign("s", csample(), eval.env = bootenv, assign.env = bootenv)
      delayedAssign("bootx", mm$x[s, , drop = FALSE], eval.env = bootenv, assign.env = bootenv)
      delayedAssign("booty", mm$y[s, , drop = FALSE], eval.env = bootenv, assign.env = bootenv)
      delayedAssign("bootivy", mm$ivy[s, , drop = FALSE], eval.env = bootenv, assign.env = bootenv)

      delayedAssign("mms",
        {
          mm1 <- list()
          mm1$extra <- mm$extra
          mm1$extra$cluster <- lapply(mm$extra$cluster, function(f) f[s])
          mm1$extra$naact <- NULL
          mm1$x <- bootx
          mm1$y <- booty
          if (!is.null(mm$ivx)) mm1$ivx <- mm$ivx[s, , drop = FALSE]
          if (!is.null(mm$ivy)) mm1$ivy <- bootivy
          if (!is.null(mm$ctrl)) mm1$ctrl <- mm$ctrl[s, , drop = FALSE]
          if (!is.null(mm$fl)) mm1$fl <- lapply(mm$fl, function(f) factor(f[s]))
          if (!is.null(weights)) mm1$weights <- mm$weights[s]
          mmdemean(mm1)
        },
        eval.env = bootenv,
        assign.env = bootenv
      )
      delayedAssign("est",
        felm.mm(mms, bootstat, exactDOF, keepX, keepCX, keepModel, kclass, fuller, onlyse, psdef = psdef),
        eval.env = bootenv,
        assign.env = bootenv
      )
      delayedAssign("beta", coef(est), assign.env = bootenv, eval.env = bootenv)
      for (n in colnames(mm$x)) {
        do.call(delayedAssign, list(n, bquote(est$coefficients[.(n), ]),
          eval.env = bootenv,
          assign.env = bootenv
        ))
        do.call(delayedAssign, list(paste0("var.", n), bquote(bootx[, .(n)]),
          eval.env = bootenv,
          assign.env = bootenv
        ))
      }

      for (n in colnames(mm$y)) {
        do.call(delayedAssign, list(n, bquote(booty[, .(n)]),
          eval.env = bootenv,
          assign.env = bootenv
        ))
      }

      if (!is.null(mm$ivy)) {
        # it's an IV-estimation, make provisions for using the OLS-version
        delayedAssign("olsmms",
          {
            if (FALSE) s <- NULL # avoid warnings about undefined s
            # make the s by evaluating mms
            #                        mms
            mm1 <- list()
            mm1$extra <- mm$extra
            mm1$extra$ivform <- NULL
            mm1$x <- cbind(mm$x[s, , drop = FALSE], mm$ivy[s, , drop = FALSE])
            mm1$y <- mm$y[s, , drop = FALSE]
            if (!is.null(mm$ctrl)) mm1$ctrl <- mm$ctrl[s, , drop = FALSE]
            if (!is.null(mm$fl)) mm1$fl <- lapply(mm$fl, function(f) factor(f[s]))
            if (!is.null(weights)) mm1$weights <- mm$weights[s]
            mmdemean(mm1)
          },
          eval.env = bootenv,
          assign.env = bootenv
        )

        delayedAssign("ols",
          felm.mm(olsmms, bootstat, exactDOF, keepX, keepCX, keepModel, onlyse = onlyse, psdef = psdef),
          eval.env = bootenv,
          assign.env = bootenv
        )

        for (n in colnames(mm$ivy)) {
          do.call(delayedAssign, list(n, bquote(bootivy[, .(n)]),
            eval.env = bootenv,
            assign.env = bootenv
          ))
        }
      }
      eval(bootexpr, bootenv)
    }
    z$boot <- replicate(Nboot, bootfun())
  }
  z
}

felm.mm <- function(mm, nostats, exactDOF, keepX, keepCX, keepModel, kclass = NULL, fuller = NULL, onlyse = FALSE, psdef = FALSE) {
  ivform <- mm$ivform
  if (is.null(ivform)) {
    # no iv, just do the thing
    z <- newols(mm, nostats = nostats[1], exactDOF = exactDOF, onlyse = onlyse, psdef = psdef)
    if (keepX) z$X <- if (is.null(mm$orig)) mm$x else mm$orig$x
    if (keepCX) {
      z$cX <- mm$x
      z$cY <- mm$y
    }
    if (keepModel) z$model <- mm$model else z$model <- mm$mfcall
    z$call <- match.call()
    return(z)
  }


  if (length(nostats) == 2) {
    nost1 <- nostats[2]
  } else {
    nost1 <- nostats[1]
  }

  ########### Instrumental variables ############


  fitnames <- makefitnames(colnames(mm$ivy))
  # should we do k-class estimation?
  if (is.null(kclass) || is.numeric(kclass)) {
    kappa <- kclass
  } else {
    KN <- ncol(mm$ivx)
    LN <- ncol(mm$x)
    N <- nrow(mm$x)
    # todo: liml

    kappa <- switch(kclass,
      `2sls` = ,
      tsls = 1.0,
      nagar = 1 + (KN - 2) / N,
      b2sls = ,
      btsls = 1 / (1 - (KN - 2) / N),
      mb2sls = ,
      mbtsls = (1 - LN / N) / (1 - KN / N - LN / N),
      liml = limlk(mm),
      fuller = limlk(mm) - fuller / (N - KN),
      stop("Unknown k-class: ", kclass, call. = FALSE)
    )
    if (identical(kclass, "liml") && fuller != 0) {
      kappa <- kappa - fuller / (N - KN)
    }
  }
  # if k-class, we should add all the exogenous variables
  # to the lhs in the 1st stage, and obtain all the residuals
  # of the instruments.  A fraction (1-kappa) of the residuals
  # are added to the fitted values when doing the 2nd stage.
  # nah, we should project on P_{Z,W}. Now, P_{Z,W} W = W
  if (!is.null(kappa)) {
    mm2 <- mm1 <- mm[names(mm) %in% c(
      "fl", "terms", "cluster", "numctrl",
      "hasicpt", "na.action", "contrasts",
      "weights"
    )]
    nmx <- colnames(mm$x)
    mm1$y <- cbind(mm$x, mm$ivy)
    mm1$x <- cbind(mm$x, mm$ivx)

    mm2$y <- mm$y
    mm2$x <- mm1$y
    mm2$orig$x <- cbind(mm$orig$x, mm$orig$ivx)
    mm2$orig$y <- cbind(mm$orig$y, mm$orig$ivy)
    #    rm(mm)
    z1 <- newols(mm1, nostats = nost1, onlyse = onlyse, psdef = psdef)

    mm2$noinst <- z1$residuals
    rm(mm1)
    #    setdimnames(mm2$x, list(NULL, c(fitnames,nmx)))
    z2 <- newols(mm2, exactDOF = exactDOF, kappa = kappa, nostats = nostats[1], onlyse = onlyse, psdef = psdef)
    if (keepX) z2$X <- if (is.null(mm2$orig)) mm2$x else mm2$orig$x
    if (keepCX) {
      z2$cX <- mm2$x
      z2$cY <- mm2$y
    }
    if (keepModel) z2$model <- mm$model else z2$model <- mm$mfcall
    z2$call <- match.call()
    return(z2)
  }


  # Now, we must build a model matrix with the endogenous variables
  # on the left hand side, the exogenous and instruments on the rhs.
  # we have already centred everything in mm. However, we must
  # rearrange it.
  # in the first stage we should have the iv left hand side on
  # the lhs, the exogenous and instruments on the rhs.
  # mm$x is an ok rhs. The y must be replaced by the ivy

  ivars <- colnames(mm$ivx)
  exlist <- colnames(mm$x)

  mm1 <- mm[names(mm) %in% c(
    "fl", "terms", "cluster", "numctrl",
    "hasicpt", "na.action", "contrasts",
    "weights"
  )]
  mm1$y <- mm$ivy
  mm1$x <- cbind(mm$x, mm$ivx)
  mm1$orig$y <- mm$orig$ivy
  mm1$orig$x <- cbind(mm$orig$x, mm$orig$ivx)
  z1 <- newols(mm1, nostats = nost1, exactDOF = exactDOF, onlyse = onlyse, psdef = psdef)
  if (keepX) z1$X <- if (is.null(mm1$orig)) mm1$x else mm1$orig$x
  if (keepCX) {
    z1$cX <- mm1$x
    z1$cY <- mm1$y
  }
  rm(mm1)

  if (!nost1) {
    z1$iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1, ivars, lhs = lh))
    names(z1$iv1fstat) <- z1$lhs
    z1$rob.iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1, ivars, type = "robust", lhs = lh))
    names(z1$rob.iv1fstat) <- z1$lhs
  }

  z1$instruments <- ivars
  z1$centred.exo <- mm$x
  z1$ivx <- mm$ivx
  z1$ivy <- mm$ivy

  # then second stage.  This is a bit more manipulation
  # We must make an mm2 with the original lhs, and with
  # the exogenous variables plus the the predicted endogenous from z1 on the rhs
  # we must set the names of the exogenous variables

  mm2 <- mm[names(mm) %in% c(
    "fl", "terms", "cluster", "numctrl", "hasicpt",
    "na.action", "contrasts", "TSS", "weights"
  )]

  mm2$x <- cbind(mm$x, z1$c.fitted.values)
  setdimnames(mm2$x, list(NULL, c(exlist, fitnames)))
  mm2$y <- mm$y
  if (!is.null(mm$orig)) {
    mm2$orig <- list(x = cbind(mm$orig$x, z1$fitted.values), y = mm$orig$y)
    setdimnames(mm2$orig$x, list(NULL, c(exlist, fitnames)))
  }

  #  rm(mm)  # save some memory

  z2 <- newols(mm2, stage1 = z1, nostats = nostats[1], exactDOF = exactDOF, kappa = kappa, onlyse = onlyse, psdef = psdef)
  if (keepX) z2$X <- if (is.null(mm2$orig)) mm2$x else mm2$orig$x
  if (keepCX) {
    z2$cX <- mm2$x
    z2$cY <- mm2$y
  }
  if (keepModel) z2$model <- mm$model else z2$model <- mm$mfcall
  rm(mm2)

  z2$call <- match.call()
  z2
}

#' Check if formula contains redundant FEs
#'
#' Print warning when there is something like fe1+fe1:fe2
#' @param formula a formula
#' @examples
#' check_redundant_fe(y ~ x)
#' check_redundant_fe(y ~ x | fe1 + fe2)
#' check_redundant_fe(y ~ x | fe1 + fe2:fe1)
#' check_redundant_fe(y ~ x | fe2 * fe1)
#' @author Grant McDermott, small adaptation by Matthieu Stigler
#' @noRd
check_redundant_fe <- function(formula) {
  fml_chk <- Formula::Formula(formula)
  has_FE <- length(attr(fml_chk, "rhs")) > 1 && !is.null(attr(fml_chk, "rhs")[[2]])
  if (has_FE) {
    fes_chk <- attr(terms(formula(fml_chk, rhs = 2)), "term.labels")
    if (any(duplicated(unlist(strsplit(fes_chk, ":"))))) {
      warning(paste(
        "Duplicated terms detected in the fixed effects slot.",
        "If you are interacting factor variables, consider excluding",
        "the parents terms, since these strictly nest the interactions",
        "and are thus redundant. See ?felm 'Details'.\n"
      ))
    }
  }
}
