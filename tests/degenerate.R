library(lfe)
set.seed(127)
options(lfe.threads=1,warn=1,digits=5)
x <- rnorm(2000,mean=2000)
x2 <- rnorm(length(x))
x3 <- 1.2*x + 0.9*x2 
x4 <- 0.8*x2 + 0.3*x3
## create individual and firm
id <- factor(sample(12,length(x),replace=TRUE))
firm <- factor(sample(7,length(x),replace=TRUE))

# these are constant on the levels
x5 <- rnorm(nlevels(id))[id]
x6 <- rnorm(nlevels(firm))[firm]
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + 0.5*x3 + 0.17*x4 + 0.8*x5 -0.3*x6 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate
est <- felm(y ~ x+x2 + x3 + x4 + x5 + x6 | id + firm)
## extract the group fixed effects
alpha <- getfe(est)
summary(est)  
#alpha
#summary(lm(y ~ x + x2 + x3 + x4 + x5 + x6 + id + firm)) # remove from cran
# merge back
ideff <- alpha[paste('id',id,sep='.'),'effect']
firmeff <- alpha[paste('firm',firm,sep='.'),'effect']

## verify that id and firm coefficients are 1
cat('accuracy:',sprintf('%.8e',coef(lm(y ~ x + x2 + x3 + x4 + x5 + x6 + ideff + firmeff-1))[7:8],'\n'))


