library(lfe)
set.seed(43)
options(lfe.threads=2,digits=5,warn=1)

g1 <- 80
g2 <- 20
g3 <- 12
N <- 1000
clu1 <- sample(g1,N, replace=TRUE)
clu2 <- (clu1 + sample(7,N,replace=TRUE)-1) %% g2
clu3 <- (clu2 + sample(3,N,replace=TRUE)-1) %% g3
clu1 <- factor(clu1)
clu2 <- factor(clu2)
clu3 <- factor(clu3)
# group specific covariate effects
ceff1 <- rnorm(nlevels(clu1), sd=0.5)[clu1]
ceff2 <- rnorm(nlevels(clu2), sd=0.4)[clu2]
ceff3 <- rnorm(nlevels(clu3))[clu3]

 # group specific errors
err1 <- rnorm(nlevels(clu1), sd=0.8)[clu1]
err2 <- rnorm(nlevels(clu2))[clu2]
err3 <- rnorm(nlevels(clu3), sd=0.5)[clu3]

x1 <- ceff1 + 0.3*ceff2 + rnorm(N)
x2 <- ceff2 + 0.2*ceff3 + rnorm(N)
x3e <- ceff3 + 0.2*(ceff2+ceff1) + rnorm(N)

f1 <- factor(sample(8,N,replace=TRUE))
x3 <- as.vector(as(f1,'sparseMatrix') %*% x3e)[f1]/tabulate(f1)[f1]
err <- err1 + err2 + err3 + abs(x1+x2*x3)*rnorm(N)
y <- x1 + x2 + x3 + err
data <- data.frame(y,x1,x2,x3,f1,clu1,clu2,clu3)
clu <- list('clu1', 'clu2', 'clu3')
summary(felm(y ~ x1 + x2 + f1|0|0|clu1+clu2+clu3, data))
#gclu <- structure(clu, method='gaure')
#summary(felm(y ~ x1 + x2 + f1|0|0|clu1+clu2+clu3, data, cmeth='gaure'))
#summary(est <- felm(y ~ x1 + x2 | f1, data, clustervar=gclu))
#ef <- structure(function(x,addnames) {
#    c(x[1],x[2:8]-x[1])
#}, verified=TRUE)
#getfe(est,ef=ef,se=TRUE, bN=200)
