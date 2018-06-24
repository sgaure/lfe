library(lfe)
options(lfe.threads=1,digits=3,warn=1)
set.seed(6553)
# single fixed effect, special case which we risk destroying when we optimize, so check it
x <- rnorm(2000)
x2 <- rnorm(length(x))
x3 <- rnorm(length(x))
## create individual and firm
id <- factor(sample(1500,length(x),replace=TRUE))
nlevels(id)
## effects
id.eff <- rnorm(nlevels(id))

## left hand side
y <- x + 0.25*x2 + 0.5*x3 + id.eff[id] + rnorm(length(x))

## estimate
est <- felm(y ~ x+x2 + x3 |id)

## extract the group fixed effects
fe <- getfe(est, se=TRUE)
## merge back
head(fe)
ideff <- fe[paste('id',id,sep='.'),'effect']

## verify that id and firm coefficients are 1
options(scipen=8)
lm(y ~ x + x2 + x3 + ideff -1)

# no factor
felm(y ~ x + x2 + x3)

# no covariate
est <- felm(y ~ 0|id)
head(getfe(est, se=TRUE))

