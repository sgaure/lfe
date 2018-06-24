library(lfe)
options(lfe.threads=2,digits=5,warn=1,lfe.eps=1e-5)
set.seed(42)
x <- rnorm(400)
x2 <- rnorm(length(x))

id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1,1)))
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
x3 <- 0.2*x + 0.3*x2 + rnorm(length(x))
x4 <- 0.1*x - 0.2*x2 + rnorm(length(x))
Q <- 0.3*x3 + 0.4*x4 + x + 0.2*x2 + id.eff[id] + 0.15*u + rnorm(length(x),sd=0.2)
R <- 0.3*x3 + 0.33*x4 + 0.2*x + 0.5*x2 + 0.7*id.eff[id] - 0.11*u + rnorm(length(x),sd=0.2)
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + Q + R + u

## estimate and print result
est <- felm(y ~ x+x2 | id+firm |(Q|R~x3+x4))
summary(est,robust=TRUE)
summary(felm(y ~ x+x2 | id+firm |(Q|R~x3+x4), kclass='liml'))
update(est, . ~ x)
# try it from within a function 
fr <- data.frame(y,x,id,firm,Q,R,x3,x4)
fun <- function() {
  Y <- y
  S <- Q
  clu <- factor(sample(10,length(x), replace=TRUE))
  felm(Y ~ x+x2 | id + firm |(Q|R ~ x3+x4), cluster=clu)
  fr <- data.frame(y,x,x2,id,firm,Q,R,x3,x4,clu)
  # test whether it finds names in the wrong place.
  `S(fit)` <- as.name('a')
  R <- as.name('b')
  felm(y ~ x+x2 | id + firm |(S|R ~ x3+x4)|clu, data=fr)
}
est <- fun()
rm(x2)
summary(est, robust=TRUE)
for(lh in est$stage1$lhs) print(summary(est$stage1, lhs=lh))
print(condfstat(est,NULL,quantiles=c(0.1,0.5,0.9)))
