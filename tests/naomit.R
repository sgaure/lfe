library(lfe)
options(lfe.threads=2,digits=5,warn=1)
set.seed(65318)
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(6,length(x),replace=TRUE,prob=c(2,rep(1,5))))

## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))
x[sample(500,20)] <- NA
y[sample(500,20)] <- NA
## estimate
summary(est <- felm(y ~ x+x2 | id+firm))
getfe(est)
summary(lm(y ~ x + x2 + id + firm -1))
