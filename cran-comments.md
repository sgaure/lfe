## Test environments
* local R installation, R 4.0.3 on Ubuntu 20.04
* win-builder: devel, release, oldrelease

## R CMD check results

0 errors | 0 warnings | 0 note

## Package adopted

* This package was adopted due to its popularity. 
* `Sprintf` issue was corrected
* Changes kindly made by CRAN maintainers were backported (commit https://github.com/MatthieuStigler/lfe/commit/dbbcb4eaea749dae2dcd071e1857308dd8e071f6)
* DOI were added in DESCRIPTION
* I cannot promise to rapidly solve the compiler flag note on Solaris. The complex code goes way beyond my understanding, but I will do my best to seek help from more qualified people. Given the popularity of the package (7000-9000 downloads per month), and given that this only triggers a note on a rather exotic platform, we were hoping you would still accept this submission? Thanks for your understanding!
