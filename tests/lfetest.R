library(lfe)
options(lfe.threads=2,digits=3,warn=1)
set.seed(655318)
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(400,length(x),replace=TRUE))
firm <- factor(sample(300,length(x),replace=TRUE,prob=c(2,rep(1,299))))

## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

# make a data frame
fr <- data.frame(y,x,x2,id,firm)
## estimate and print result
print(est <- felm(y ~ x+x2|id+firm, data=fr))

## extract the group fixed effects
tail(getfe(est))

head(model.matrix(est))
head(model.matrix(est,centred=FALSE))
head(model.matrix(felm(y ~ x+x2|id+firm, data=fr,keepCX=TRUE)))
head(model.matrix(felm(y ~ x+x2|id+firm, data=fr,keepX=TRUE), centred=FALSE))
