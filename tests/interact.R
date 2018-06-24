library(lfe)
options(lfe.threads=2,digits=3,warn=1)
set.seed(42)
x <- rnorm(5000,mean=2)
x2 <- rnorm(length(x))
x3 <- rexp(length(x))
#ix <- rep(2.0,length(x))
## create individual and firm
id <- factor(sample(150,length(x),replace=TRUE))
ix <- matrix(rnorm(2*length(x)),length(x),2)*sqrt(as.numeric(id)/nlevels(id))
colnames(ix) <- c('col1','col2')
firm <- factor(sample(130,length(x),replace=TRUE))
fx <- rnorm(length(x))*sqrt(as.numeric(id)/nlevels(firm))
shoe <- factor(sample(10,length(x),replace=TRUE))
shirt <- factor(sample(10,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
shoe.eff <- rnorm(nlevels(shoe))
shirt.eff <- rnorm(nlevels(shirt))
## left hand side
y <- x + 0.25*x2 + 0.5*x3 + rowSums(ix*id.eff[id]) + id.eff[id] + fx*firm.eff[firm] + firm.eff[firm] + shoe.eff[shoe] + shirt.eff[shirt] + rnorm(length(x))

## estimate
summary(est <- felm(y ~ x+x2+x3 | ix:id+id+fx:firm+firm+shoe+shirt, exactDOF=T))
getfe(est)

# make sure this one works, 
summary(est <- felm(y ~ x | x2:id + x3:id + id))


# compare with lm
x <- rnorm(100)
id <- factor(sample(5,length(x), replace=T))
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
firm <- factor(sample(7,length(x), replace=T))
fx <- rnorm(length(x))
ix <- matrix(rnorm(2*length(x)),length(x),2)
y <- x + rowSums(ix*id.eff[id]) + fx*firm.eff[firm] + firm.eff[firm] + rnorm(length(x))
summary(est <- felm(y ~ x | ix:id + fx:firm + firm, exactDOF='rM'))
getfe(est)
summary(lm(y ~ x + ix:id + fx:firm + firm))

