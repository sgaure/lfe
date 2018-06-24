library(lfe)
# From http://diffuseprior.wordpress.com/2012/06/15/standard-robust-and-clustered-standard-errors-computed-in-r/
set.seed(123)
options(lfe.threads=2,digits=5,warn=1)
ols <- function(form, data, robust=FALSE, cluster=NULL,digits=getOption('digits')){
    r1 <- lm(form, data)
      if(length(cluster)!=0){
            data <- na.omit(data[,c(colnames(r1$model),cluster)])
                r1 <- lm(form, data)
          }
      X <- model.matrix(r1)
      n <- dim(X)[1]
      k <- dim(X)[2]
      if(robust==FALSE & length(cluster)==0){
            se <- sqrt(diag(solve(crossprod(X)) * as.numeric(crossprod(resid(r1))/(n-k))))
                res <- cbind(coef(r1),se)
          }
      if(robust==TRUE){
            u <- matrix(resid(r1))
                meat1 <- t(X) %*% diag(diag(crossprod(t(u)))) %*% X
                dfc <- n/(n-k)
                se <- sqrt(dfc*diag(solve(crossprod(X)) %*% meat1 %*% solve(crossprod(X))))
                res <- cbind(coef(r1),se)
          }
      if(length(cluster)!=0){
            clus <- cbind(X,data[,cluster],resid(r1))
                colnames(clus)[(dim(clus)[2]-1):dim(clus)[2]] <- c(cluster,"resid")
                m <- dim(table(clus[,cluster]))
                dfc <- (m/(m-1))*((n-1)/(n-k))
                uclust  <- apply(resid(r1)*X,2, function(x) tapply(x, clus[,cluster], sum))
                se <- sqrt(diag(solve(crossprod(X)) %*% (t(uclust) %*% uclust) %*% solve(crossprod(X)))*dfc)
                res <- cbind(coef(r1),se)
          }
      res <- cbind(res,res[,1]/res[,2],(1-pnorm(abs(res[,1]/res[,2])))*2)
      res1 <- matrix(as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),res)),nrow=dim(res)[1])
      rownames(res1) <- rownames(res)
      colnames(res1) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
      return(res1)
  }



x <- rnorm(1000) 
f1 <- sample(8,length(x), repl=T)
clu <- factor(sample(10,length(x), replace=T))
cluerr <- rnorm(nlevels(clu))[clu]
clu2 <- factor(sample(10,length(x), replace=T))
cluerr2 <- rnorm(nlevels(clu2))[clu2]
err <- abs(x)*rnorm(length(x)) + cluerr + cluerr2
y <- x +rnorm(nlevels(clu),sd=0.3)[clu] +  log(f1) + err
dat <- data.frame(y, x, f1=factor(f1), cluster=clu,cluster2=clu2)
summary(felm(y ~x |f1, dat))
# CGM clustering, i.e. one factor means standard one-way clustering
summary(felm(y ~x + f1, dat, clustervar='clu'))
# this will make my experimental clustered errors for f1, typically better for few groups
# summary(felm(y ~x + f1|0|0|cluster+cluster2, dat))
summary(felm(y ~x + f1|0|0|cluster+cluster2, dat, psdef=FALSE))
# this will sample them for f1, also test having cluster in the third component
summary(estg <- felm(y ~x | f1|0|cluster, dat))
# Comparable estimable function
ef <- function(gamma, addnames) {
  ref1 <- gamma[[1]]
  res <- c(gamma[[1]],gamma[2:8]-gamma[[1]])
  if(addnames) {
    names(res) <- c('icpt',paste('f1',2:8,sep='.'))
  }
  res
}
getfe(estg,ef=ef,se=TRUE,bN=200)

#summary(estr <- felm(y ~x + G(f1) + G(f2), dat), robust=TRUE)
#ols(y ~x + f1 + f2, dat, robust=TRUE)
#getfe(estr,ef=ef,se=T,bN=2000, robust=TRUE)
#ols(y ~x + f1 + f2, dat, cluster="cluster")
