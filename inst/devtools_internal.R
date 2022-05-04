

## Initiate some docs
# usethis::use_cran_comments(open = FALSE)

## Spell check
devtools::spell_check()



### Tests
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()


## Check on R-hub
plats_default <- c("windows-x86_64-devel", "ubuntu-gcc-release", "fedora-clang-devel", 
                   "linux-x86_64-rocker-gcc-san") # from: rhub:::default_cran_check_platforms(".")
plats_add <- c("solaris-x86-patched", "solaris-x86-patched-ods")
plats_all <-  c(plats_default, plats_add)
devtools::check_rhub(platform = plats_all,
                     env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false"), ##testhat not working on Solaris by May 2022
                     interactive = FALSE)

## reverse dependencies
# devtools::install_github("r-lib/revdepcheck")
revdepcheck::revdep_check(num_workers = 1)

### Upload to CRAN
devtools::release_checks()
devtools::release()
# direct: devtools::submit_cran()