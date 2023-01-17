library(lfe)
options(lfe.threads=2,digits=5,warn=1)
set.seed(655320)
x <- rnorm(5000,mean=200)
x2 <- rnorm(length(x))
x3 <- rexp(length(x))
## create individual and firm
id <- factor(sample(1500,length(x),replace=TRUE))
firm <- factor(sample(1300,length(x),replace=TRUE))
shoe <- factor(sample(100,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
shoe.eff <- rnorm(nlevels(shoe))
## left hand side
y <- x + 0.25*x2 + 0.5*x3 + id.eff[id] + firm.eff[firm] + shoe.eff[shoe] + rnorm(length(x))

## estimate
summary(est <- felm(y ~ x+x2 + x3 | id + firm + shoe))
cat('Components:',nlevels(est$cfactor),'largest:',sum(est$cfactor == '1'),'\n')
## extract the group fixed effects
for(ef in c('ln','ref')) {
  fe <- getfe(est,ef=ef)
  ## merge back

  ideff <- fe[paste('id',id,sep='.'),'effect']
  firmeff <- fe[paste('firm',firm,sep='.'),'effect']
  shoeeff <- fe[paste('shoe',shoe,sep='.'),'effect']

  ## verify that id and firm coefficients are 1
  options(scipen=8)
  print(summary(lm(y ~ x + x2 + x3 + ideff + firmeff + shoeeff -1),digits=8))
}

# Perform a bootstrap
a <- felm(y ~ x+x2 + x3 | id + firm, nostats=TRUE, Nboot=10, bootexpr=quote(x/x3*x2))
mean(a$boot)
sd(a$boot)
