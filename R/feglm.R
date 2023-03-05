#' @title Fixed-Effects Poisson Pseudo Maximum Likelihood (PPML)
#' @inheritParams felm
#' @param offset EXPLAIN
#' @param tol EXPLAIN
#' @importFrom Formula Formula
#' @importFrom Matrix Diagonal
#' @export fepois
fepois <- function(formula, data,
                         offset = NULL,
                         subset = NULL,
                         robust = TRUE,
                         cluster = NULL,
                         pseudo_r2 = TRUE,
                         tol = 1e-10) {
  if (!is.null(subset)) { data <- data[subset, ] }
  
  if (is.character(formula)) { formula <- as.formula(formula) }
  if (is.character(cluster)) { cluster <- as.formula(paste0("~", cluster)) }
  
  formula <- Formula(formula)
  
  offset2 <- offset
  if (is.null(offset)) offset <- rep(0, nrow(data))
  
  vardep <- all.vars(formula(formula, lhs = 1, rhs = 0))
  vardep <- data[[vardep]]
  
  if (min(vardep) < 0) {
    stop("y should be greater or equals to zero.")
  }
  
  fe <- all.vars(formula(formula, lhs = 0, rhs = 2))
  
  for (f in fe) {
    if (!is(data[, f], "factor")) {
      data[, f] <- as.factor(data[, f])
    }
  }
  
  max_vardep <- max(vardep)
  vardep <- vardep / max_vardep
  
  mu <- (vardep + 0.5) / 2
  eta <- log(mu) - offset
  z <- eta + (vardep - mu) / mu
  
  # Formula

  varind <- all.vars(formula(formula, lhs = 0, rhs = 1))
  
  formula <- as.formula(paste0(
    "z ~ ", ifelse(is.null(varind), " -1 ", varind), " | ",
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
        "| ", fe_tmp, " | 0 | 0"
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
    
    cluster_name <- cluster
    
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
    x_var <- as.matrix(data[, varind])
    beta <- as.matrix(reg$coefficients)
    x_beta <- exp(x_var %*% reg$coefficients + offset + x_fe)
  } else {
    x_beta <- exp(offset + x_fe)
  }
  
  reg$fitted.values <- x_beta
  
  if (isTRUE(pseudo_r2)) {
    reg$pseudo_rsq <- cor(data[, vardep, drop = TRUE], x_beta[, "z", drop = T])^2
  }
  
  class(reg) <- "gravity.ppml"
  return(reg)
}

# summary.gravity.ppml <- function(object) {
#   class(object) <- "summary.gravity.ppml"
#   object
# }
# 
# print.summary.gravity.ppml <- function(object) {
#   cat("Coefficients: \n")
#   results <- data.frame(
#     Estimate = object$coefficients,
#     `Std. Error` = object$se,
#     `t value` = object$tval,
#     `Pr(>|t|)` = object$pval
#   )
#   results <- as.matrix(results)
#   colnames(results) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#   printCoefmat(results, digits = 4)
# }
# 
# predict.gravity.ppml <- function(object, newdata = NULL, offset = NULL) {
#   if (is.null(offset)) offset <- rep(0, nrow(newdata))
#   
#   fe <- names(object$fe)
#   x_fe <- newdata[, fe]
#   x_fe$order <- 1:nrow(x_fe)
#   len_fe <- length(fe)
#   
#   for (i in 1:len_fe) {
#     fe_tmp <- getfe(object)
#     fe_tmp <- fe_tmp[fe_tmp$fe == fe[i], c("idx", "effect")]
#     
#     colnames(fe_tmp) <- c(
#       fe[i],
#       paste0("fe_", fe[i])
#     )
#     
#     x_fe <- merge(x_fe,
#                              fe_tmp,
#                              by = fe[i],
#                              all.x = TRUE
#     )
#   }
#   
#   x_fe <- x_fe[
#     order(x_fe$order),
#     -(len_fe + 1)
#   ]
#   
#   x_fe[, 1:len_fe] <- sapply(x_fe[, 1:len_fe], as.character)
#   object$fixed.effects <- x_fe
#   
#   x_fe <- x_fe[, !names(x_fe) %in% fe]
#   x_fe <- apply(x_fe, 1, sum)
#   
#   x <- rownames(object$beta)
#   
#   if (!is.null(x)) {
#     x_var <- as.matrix(newdata[, x])
#     beta <- as.matrix(object$coefficients)
#     x_beta <- exp(x_var %*% object$coefficients + offset + x_fe)
#   } else {
#     x_beta <- exp(offset + x_fe)
#   }
#   
#   x_beta
# }
