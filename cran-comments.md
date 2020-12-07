## Test environments
* local R installation, R 4.0.3 on Ubuntu 20.04
* win-builder: devel, release, oldrelease

## R CMD check results

0 errors | 0 warnings | 0 note

## Package adopted

* This package was temporarily adopted due to its popularity. 
* I was unfortunately not able to solve the issue `Compilation used the following non-portable flag(s): ‘-mt’`. Compiling code in lfe is very complex (see `configure` and `tools/ax_pthread`), and goes beyond my understanding and that of multiple users (see Github discussion where I asked for help: https://github.com/sgaure/lfe/issues/41#issuecomment-739369852). Given the popularity of the package (7000-9000 downloads per month), and given that this only triggers a note on a rather exotic platform, we were hoping you would still accept this submission? Thanks for your understanding!
