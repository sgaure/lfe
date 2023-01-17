
context("Test clustering comparisons of https://github.com/sgaure/lfe/issues/1")

test_that("two-level clustering works", {
  # Does this regression make sense? Definitely not!
  cgm_est <- stats::coef(summary(felm(
    Ozone ~ Solar.R + Wind | Day | 0 | Day + Month, data = airquality
  )))
  # Results from running the above with
  # https://github.com/grantmcdermott/lfe/tree/66f9c0fab184b366ed60b9732475b58a6acee5dd
  cgm_expected <- matrix(c(
    #  Coef        Cluster SE   t-stat    p-value
       0.09926081, 0.02835782,  3.500297, 0.024889465,
      -5.57693782, 1.03966448, -5.364171, 0.005830014
    ),
    byrow=TRUE, nrow=2,
    dimnames=list(
      c("Solar.R", "Wind"), # rows
      c("Estimate", "Cluster s.e.", "t value", "Pr(>|t|)") # columns
    )
  )
  expect_equal(cgm_est, cgm_expected, tolerance = 0.0001, scale = 1)

  reghdfe_stats <- stats::coef(summary(felm(
    Ozone ~ Solar.R + Wind | Day | 0 | Day + Month,
    data = airquality,
    cmethod = "reghdfe"
  )))
  # Results from running the above in Stata's reghdfe
  # To export:
  # names(airquality) <- gsub('.', '_', names(airquality), fixed=TRUE)
  # haven::write_dta(data=airquality, path="airquality.dta", version=14)
  # In stata:
  # use airquality, clear
  # reghdfe Ozone Solar_R Wind , absorb(Day) cluster(Month Day)
  # regsave, tstat pval
  # list
  #    +-----------------------------------------------------------------------+
  #    |     var        coef     stderr      tstat       pval     N         r2 |
  #    |-----------------------------------------------------------------------|
  # 1. | Solar_R    .0992608   .0285888    3.47202   .0255366   109   .6644437 |
  # 2. |    Wind   -5.576938   1.053322   -5.29462   .0061094   109   .6644437 |
  # 3. |   _cons     78.6738   10.32352   7.620831   .0015917   109   .6644437 |
  #    +-----------------------------------------------------------------------+
  # Results above are with Stata 15.1 and reghdfe 5.6.8 03mar2018
  reghdfe_expected <- matrix(c(
    #  Coef       Cluster SE  t-stat   p-value
       0.0992608, 0.0285888,  3.47202, 0.0255366,
      -5.576938,  1.053322,  -5.29462, 0.0061094
    ),
    byrow=TRUE, nrow=2,
    dimnames=list(
      c("Solar.R", "Wind"), # rows
      c("Estimate", "Cluster s.e.", "t value", "Pr(>|t|)") # columns
    )
  )
  expect_equal(reghdfe_stats, reghdfe_expected,tolerance = 0.001, scale = 1)
})
