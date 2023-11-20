file tinydata.csv
vars x x2 year id firm y ife ffe yfe
model y ~ x + x2 + year + G(id)+G(firm)
dummy year
#tol 1e-3
# nofe
merge





