

## Initiate some docs
# usethis::use_news_md()
# usethis::use_cran_comments(open = FALSE)

### Tests
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

devtools::check_rhub(interactive = FALSE)

### Upload to CRAN
devtools::release_checks()
devtools::release()