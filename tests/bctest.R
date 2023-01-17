library(lfe)
options(lfe.threads=1,digits=4,warn=1)
set.seed(42)
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(40,length(x),replace=TRUE))
firm <- factor(sample(30,length(x),replace=TRUE,prob=c(2,rep(1,29))))
foo <- factor(sample(20,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
foo.eff <- rnorm(nlevels(foo))
## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + foo.eff[foo] + rnorm(length(x))

# make a data frame
fr <- data.frame(y,x,x2,id,firm,foo)
## estimate and print result
est <- felm(y ~ x+x2|id+firm+foo, data=fr, keepX=TRUE)

alpha=getfe(est)
bccorr(est,alpha,corrfactors=c(3,1))
fevcov(est,alpha)
varvars(est,alpha)

