#' @title Fixed-Effects Poisson Pseudo Maximum Likelihood (PPML)
#' @inheritParams felm
#' @param offset this can be used to specify an \emph{a priori} known component
#'  to be included in the linear predictor during fitting. This should be
#'  \code{NULL} or a numeric vector or matrix of extents matching those of the
#'  response. One or more \code{\link{offset}} terms can be included in the
#'  formula instead or as well, and if more than one are specified their sum is
#'  used. See \code{\link{model.offset}}.
#' @param robust logical value to return a robust standard error computation.
#' @param cluster optional variable to group by and compute sandwich-type
#' robust standard errors. Should be a formula of the form `~x_j` or
#' an object that be coerced to a formula.
#' @param tol tolerance value for GLM convergence criteria.
#' @importFrom Formula Formula
#' @importFrom Matrix Diagonal
#' @importFrom methods is
#' @seealso felm
#' @export fepois
fepois <- function(formula, data,
                         offset = NULL,
                         subset = NULL,
                         robust = TRUE,
                         cluster = NULL,
                         tol = 1e-10) {
  if (!is.null(subset)) { data <- data[subset, ] }
  
  if (is.character(formula)) { formula <- as.formula(formula) }
  if (is.character(cluster)) { cluster <- as.formula(paste0("~", cluster)) }
  if (is.character(offset)) { offset <- as.formula(paste0("~", offset)) }
  
  formula <- Formula(formula)
  
  offset2 <- offset
  
  if (is.null(offset)) {
    offset <- rep(0, nrow(data))
  } else {
    offset <- data[[all.vars(formula(offset, lhs = 1, rhs = 0))]]
  }
  
  vardep <- all.vars(formula(formula, lhs = 1, rhs = 0))
  vardep2 <- vardep # for pseudo rsq later
  vardep <- data[, vardep, drop = TRUE]
  
  if (min(vardep) < 0) {
    stop("y should be greater or equals to zero.")
  }
  
  fe <- all.vars(formula(formula, lhs = 0, rhs = 2))
  
  for (f in fe) {
    if (!is(data[, f], "factor")) {
      data[, f] <- as.factor(data[, f])
    }
  }
  
  max_vardep <- max(vardep, na.rm = TRUE)
  vardep <- vardep / max_vardep
  
  mu <- (vardep + 0.5) / 2
  eta <- log(mu) - offset
  z <- eta + (vardep - mu) / mu
  
  # Formula

  varind <- all.vars(formula(formula, lhs = 0, rhs = 1))
  
  formula <- as.formula(paste0(
    "z ~ ", ifelse(is.null(varind), " -1 ", paste0(varind, collapse = " + ")), " | ",
    paste0(fe, collapse = " + ")
  ))
  
  dif <- 1
  rss1 <- 1
  
  while (abs(dif) > tol) {
    reg <- felm(formula = formula, data = data, weights = mu)
    
    eta <- z - reg$residuals + offset
    mu <- exp(eta)
    z <- (eta - offset) + (vardep - mu) / mu
    
    res <- vardep - mu
    rss2 <- sum(res^2)
    dif <- rss2 - rss1
    rss1 <- rss2
    dev <- 2 * max_vardep * sum(vardep[vardep > 0] * log(vardep[vardep > 0] / mu[vardep > 0]))
  }
  
  z <- z + log(max_vardep)
  reg <- felm(formula = formula, data = data, weights = mu)
  
  if (!is.null(varind)) {
    z <- data.frame(id = seq_len(nrow(data)))
    for (i in varind) {
      fe_tmp <- paste0(fe, collapse = " + ")
      formula_tmp <- as.formula(paste0(
        i, " ~ -1 ",
        ifelse(!is.null(offset2), " + offset ", ""),
        "| ", fe_tmp
      ))
      fit.tmp <- felm(formula = formula_tmp, data = data, weights = mu)
      z[[i]] <- fit.tmp$residuals
    }
    z <- z[, -1]
    z <- as.matrix(z)
    
    n <- reg$N
    k <- length(varind)
    W1 <- Diagonal(mu, n = n)
    bread <- solve(t(z) %*% W1 %*% z)
    
    res <- vardep - mu
    
    if (robust) {
      if (is.null(cluster)) {
        W2 <- Diagonal((res^2), n = n)
        meat <- t(z) %*% W2 %*% z
      } else {
        cluster <- data[[all.vars(formula(cluster, lhs = 1, rhs = 0))]]
        m <- length(unique(cluster))
        dfc <- (m / (m - 1)) * ((n - 1) / (n - k))
        
        meat <- matrix(0, nrow = length(varind), ncol = length(varind))
        
        for (i in unique(cluster)) {
          z_tmp <- as.matrix(z[cluster == i, , drop = FALSE])
          res_tmp <- res[cluster == i]
          W2_tmp <- res_tmp %*% t(res_tmp)
          meat <- meat + (t(z_tmp) %*% W2_tmp %*% z_tmp)
        }
      }
    }
    
    vcov <- if (robust) {
      bread %*% meat %*% bread
    } else {
      bread
    }
    
    reg$vcv <- vcov
    reg$se <- sqrt(diag(reg$vcv))
    reg$tval <- reg$coefficients / reg$se
    reg$pval <- 1 - pnorm(abs(reg$tval))
    
    if (robust) {
      reg$vcv <- vcov
      reg$se <- sqrt(diag(reg$vcv))
      reg$tval <- reg$coefficients / reg$se
      reg$pval <- 1 - pnorm(abs(reg$tval))
      
      reg$robustvcv <- vcov
      reg$rse <- sqrt(diag(reg$robustvcv))
      reg$rtval <- reg$coefficients / reg$rse
      reg$rpval <- 1 - pnorm(abs(reg$rtval))
    }
  }
  
  # use drop = F to ensure the data frame is not converted to vector
  # https://stackoverflow.com/q/75639224/3720258
  x_fe <- data[, fe, drop = FALSE]
  len_fe <- length(fe)
  
  for (i in seq_len(len_fe)) {
    fe_tmp <- getfe(reg)
    fe_tmp <- fe_tmp[fe_tmp$fe == fe[i], c("idx", "effect")]
    
    colnames(fe_tmp) <- c(fe[i], paste0("fe_", fe[i]))
    x_fe <- merge(x_fe, fe_tmp, by = fe[i], all.x = TRUE)
  }
  x_fe[, seq_len(len_fe)] <- sapply(x_fe[, seq_len(len_fe)], as.character)
  reg$fixed.effects <- x_fe
  
  x_fe <- x_fe[, !names(x_fe) %in% fe, drop = FALSE]
  x_fe <- apply(x_fe, 1, sum)
  
  if (!is.null(varind)) {
    reg$fitted.values <- exp(as.matrix(data[, varind]) %*% reg$coefficients + offset + x_fe)
  } else {
    reg$fitted.values <- exp(offset + x_fe)
  }
  
  class(reg) <- "fepois"
  return(reg)
}

#' @exportS3Method
summary.fepois <- function(object) {
  class(object) <- "summary.fepois"
  return(object)
}

#' @exportS3Method
print.summary.fepois <- function(object) {
  cat("Coefficients: \n")
  results <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `t value` = object$tval,
    `Pr(>|t|)` = object$pval
  )
  results <- as.matrix(results)
  colnames(results) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  return(printCoefmat(results, digits = 4))
}

#' @exportS3Method
predict.fepois <- function(object, newdata = NULL, offset = NULL, type = "link") {
  stopifnot(any(type %in% c("link", "response")))
  
  if (is.null(offset)) offset <- rep(0, nrow(newdata))

  fe <- names(object$fe)
  x_fe <- newdata[, fe]
  x_fe$order <- 1:nrow(x_fe)
  len_fe <- length(fe)

  for (i in 1:len_fe) {
    fe_tmp <- getfe(object)
    fe_tmp <- fe_tmp[fe_tmp$fe == fe[i], c("idx", "effect")]

    colnames(fe_tmp) <- c(fe[i], paste0("fe_", fe[i]))

    x_fe <- merge(x_fe, fe_tmp, by = fe[i], all.x = TRUE)
  }

  x_fe <- x_fe[order(x_fe$order), -(len_fe + 1)]

  x_fe[, 1:len_fe] <- sapply(x_fe[, 1:len_fe], as.character)
  object$fixed.effects <- x_fe

  x_fe <- x_fe[, !names(x_fe) %in% fe]
  x_fe <- apply(x_fe, 1, sum)

  x <- rownames(object$beta)

  # return x_beta
  if (!is.null(x)) {
    if (type == "link") {
      return(as.matrix(newdata[, x]) %*% object$coefficients + offset + x_fe) 
    }
    if (type == "response") {
      return(exp(as.matrix(newdata[, x]) %*% object$coefficients + offset + x_fe))
    }
  } else {
    if (type == "link") {
      return(offset + x_fe)
    }
    if (type == "response") {
      return(exp(offset + x_fe)) 
    }
  }
}
