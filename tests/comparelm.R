library(lfe)
options(lfe.threads=2,digits=5,warn=1)
set.seed(65318)
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(6,length(x),replace=TRUE,prob=c(2,rep(1,5))))

## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))
y2 <- 0.1*x + 0.2*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate
smry_est <- summary(est <- felm(y ~ x+x2 | id + firm))
getfe(est)
smry_est_lm <- summary(est_lm <- lm(y ~ x + x2 + id + firm))
smry_est_lm

all.equal(coef(smry_est), coef(smry_est_lm)[c("x", "x2"),])

## vcov, confint
all.equal(vcov(est), vcov(est_lm)[c("x", "x2"), c("x", "x2")])
all.equal(confint(est), confint(est_lm, parm = c("x", "x2")))
all.equal(stats::confint.default(est), stats::confint.default(est_lm, parm = c("x", "x2")))

## should differ:
isTRUE(all.equal(confint(est, type ="iid"), confint(est, type ="robust")))

## multi response
est_felm_multi <- felm(y + y2 ~ x+x2 | id + firm)
est_lm_multi <- lm(cbind(y, y2) ~ x+x2 + id + firm)

all.equal(coef(est_felm_multi), coef(est_lm_multi)[c("x", "x2"),])
all.equal(coef(summary(est_felm_multi, lhs = "y")),
          coef(summary(est_lm_multi))[[1]][c("x", "x2"),])
all.equal(coef(summary(est_felm_multi, lhs = "y2")),
          coef(summary(est_lm_multi))[[2]][c("x", "x2"),])

all.equal(confint(est_felm_multi),
          confint(est_lm_multi)[c("y:x", "y:x2", "y2:x", "y2:x2"),])

