library(lfe)
options(lfe.threads=1,digits=4)
set.seed(42)
N <- 1000
K <- 40
id <- factor(sample(K,N, replace=TRUE))
firm <- factor(sample(K,N, replace=TRUE))
ideff <- rnorm(nlevels(id))
firmeff <- rnorm(nlevels(firm))
x <- ideff[id] + rnorm(N)
y <- x  + ideff[id] - firmeff[firm] + rnorm(N,sd=2)
w <- abs(2*x + ideff[id] + firmeff[firm]) + runif(N,0.2,0.4)

# lm says:
summary(lm(y ~ x + id + firm, wei=w))
# felm w/o projection says
summary(felm(y ~ x + id + firm, wei=w))
# felm with projection says
est <- felm(y ~x | id + firm, wei=w, keepX=TRUE)
summary(est)
head(sandwich::estfun(est))

# test getfe:
Kid <- nlevels(id)
sid <- 1:Kid
Kfirm <- nlevels(firm)
sfirm <- 1:Kfirm + Kid
ef <- function(x, addnames=FALSE) {
  icpt <- x[1] + x[Kid+1]
  x[sid[-1]] <- x[sid[-1]] - x[1]
  x[1] <- icpt
  x[sfirm] <- x[sfirm] - x[Kid+1]
  if(addnames) names(x) <- c('(Intercept)',paste('id',sid[-1],sep=''),paste('firm',1:Kid,sep=''))
  x[c(sid,sfirm[-1])]
}
head(getfe(est,ef=ef,se=TRUE))
# test fevcov
print(fv <- fevcov(est))
message('correlation:', round(cov2cor(fv)[1,2],4))
v1 <- lfe:::wvar(ideff[id],w)
v2 <- lfe:::wvar(firmeff[firm],w)
cv <- lfe:::wcov(ideff[id], -firmeff[firm], w)
message(
    'idvar ',round(v1,4),
    ' firmvar ',round(v2,4),
    ' cov ',round(cv,4),
    ' cor ',round(cv/sqrt(v1*v2),4)
    )
bccorr(est)

varvars(est)
varvars(est,biascorrect=TRUE)
