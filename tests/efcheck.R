library(lfe)
set.seed(65513)
options(lfe.threads=1,digits=4,warn=1)
x <- rnorm(52)
x2 <- rnorm(length(x))
x3 <- 0.2*x + 0.1*x2
## create individual and firm
id <- factor(sample(30,length(x),prob=c(2,2,rep(1,28)),replace=TRUE))
firm <- factor(sample(32,length(x),prob=c(2,2,rep(1,30)),replace=TRUE))
year <- factor(sample(3,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
year.eff <- rnorm(nlevels(year))
## left hand side
y <- x + 0.25*x2 + x3 + id.eff[id] + firm.eff[firm] + year.eff[year] + rnorm(length(x))

## estimate and print result
est <- felm(y ~ x+x2+x3 |id+firm+year, exactDOF=TRUE)

## extract the group fixed effects
head(getfe(est))
tail(alpha <- getfe(est,ef='ln'))
# get the names to use below, just to make it easier
# lower precision in output

nm <- rownames(alpha)
getfe(est,ef='zm',se=TRUE)
getfe(est,ef='zm2',se=TRUE)
ef <- function(v,addnames) {
  names(v) <- nm
  w <- c(v['id.2']-v['id.1'],exp(v['id.2']-v['id.1']),
         v['id.2']+v['firm.2']+v['year.1'],v['id.5']+v['firm.4']+v['year.2'])
  if(addnames) names(w) <-c('id2-id1','exp(id2-id1)','id2+f2+y1','id5+f4+y2')
  w
}
getfe(est,ef=ef,se=TRUE)

# test whether we have estimable functions

R <- est$r.residuals - est$residuals

cat('myef :',is.estimable(ef,est$fe,R),'\n')
for(n in c('ref','zm','zm2','ln')) {
  cat(n,':',is.estimable(efactory(est,n),est$fe,R),'\n')
}
