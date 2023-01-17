library(lfe)
options(digits=3)
set.seed(42)
B <- matrix(rnorm(40000),200)
A <- band(crossprod(B), -30,30)
tx <- rnorm(ncol(A))
b <- A %*% tx
fun <- function(x) A %*% x
sol <- cgsolve(fun,b,eps=-0.001, symmtest=TRUE)
sqrt(sum((sol-tx)^2))

A[5,6] <- 0.24
options(warn=2)
try(cgsolve(fun,b,symmtest=TRUE))
