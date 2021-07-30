## Test environments
* local R installation, R 4.1.0 on Ubuntu 20.04
* win-builder: devel, release, oldrelease

## R CMD check results

0 errors | 0 warnings | 0 note

## Package adopted

* This version solves the issue with the latest changes in the API with SET_REFCNT
* This version changes also a NA to NaN change in R, hence version R < 4.1 will show a warning
* There is still a compiler flag note on Solaris. Note that this package was adopted due to its popularity, as it had been archived. The complex code goes way beyond my understanding, so I cannot address the Solaris note. I hope that you will still allow to publish the package on CRAN, as it is a very popular one, and the Solaris issue is a minor issue on a rather exotic platform. Thanks!
