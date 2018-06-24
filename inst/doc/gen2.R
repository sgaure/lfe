#require('mvtnorm')
# individuals
Nind <- 15000
# firms
Nfirms <- 5000
# avg observations per individual
obs <- 15
years <- 1990:(1990+obs-1)
Nrows <- Nind*obs





#varcorr matrix of observed covariates
vcv <- matrix(c(1,0.5,0.5,0.71),ncol=2)
X <- rmvnorm(Nrows,c(0,0),vcv)
x <- X[,1]
x2 <- X[,2]

# individual effects
ife <- rchisq(Nind,2)
ife <- runif(Nind)

# normalize variance
ife <- ife/sqrt(var(ife))

# firm fixed effects
ffe <- rexp(Nfirms) 
#ffe <- runif(Nfirms)
ffe <- ffe/sqrt(var(ffe))

yfe <- seq(0.0,(obs-1)/10,length.out=obs)
#yfe[[4]] <- yfe[[4]]*100
#yfe <- rnorm(obs)


# hmm, now generate data
# start out with a random assignment to firms
# firms has a probability of being moved to
fprob <- rchisq(Nfirms,10)
#fprob <- rep(1,Nfirms)
# if want to model assortative matching, we should
# let the probability vary with the fixed effect of
# the individual and the firm.  Hmm, how to do this in a 
# simple fashion?

year = factor(rep(years,Nind))
yint <- as.integer(levels(year))[year]
tab <- data.frame(x=x,x2=x2,year=yint,id=rep(1:Nind,each=obs),firm=0,y=0)
firmy <- sample(1:Nfirms,Nind,replace=T,prob=fprob)
for(i in years) {
   cat('Doing year',i,'\n')
   tab[year == i,'firm'] <- firmy
   # Now, figure out who changes  10% prob
   ch <- runif(Nind) < 0.1
   nch <- sum(ch)
   firmy[ch] <- sample(1:Nfirms,nch,replace=T,fprob)
}

firm <- as.factor(tab[,'firm'])
id <- as.factor(tab[,'id'])

tab[,'ife'] <- ife[id]
tab[,'ffe'] <- ffe[firm]   
tab[,'yfe'] <- yfe[year]
# introduce some correlation
#tab[,'x'] <- tab[,'x'] + 0.3*ife[id] + 0.2*ffe[firm]
#tab[,'x2'] <- tab[,'x2'] + 0.2*ffe[firm] - 0.3*ife[id]
tab[,'y'] <- tab[,'ife'] + tab[,'ffe'] + tab[,'yfe'] + 0.5*x + 0.25*x2 + rnorm(Nrows,sd=sqrt(4))
write.table(tab,'tinydata.csv',row.names=FALSE,col.names=FALSE)
if(!interactive()) quit('n')
