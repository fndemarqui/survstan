

# https://community.rstudio.com/t/designing-hex-logos/2756
# https://pkgdown.r-lib.org/reference/build_favicons.html

# pacote para criar log de pacotes:
#install.packages("hexSticker")

usethis::use_git_ignore("inst/stuff")
usethis::use_build_ignore("inst/stuff")
usethis::use_git_ignore(".Rd2pdf10127")
usethis::use_build_ignore(".Rd2pdf10127")


# One time set-up
# usethis::use_mit_license()
# usethis::use_readme_rmd()
# usethis::use_news_md()
# usethis::use_cran_comments()

# One time set-up with Git, Github, and Github Actions
# gitcreds::gitcreds_set()
# usethis::use_git(message = "Initial commit")
# usethis::use_github(private = FALSE)
# usethis::use_github_action("pkgdown")
# usethis::use_pkgdown_github_pages()
# usethis::use_github_action("check-standard")



usethis::use_git_ignore("cran-comments.md")
usethis::use_build_ignore("cran-comments.md")


library(rstantools)
example(source) # defines the sourceDir() function
try(roxygen2::roxygenize(load_code = sourceDir), silent = TRUE)
# rstan_config()
pkgbuild::compile_dll(force = TRUE)

roxygen2::roxygenize()
devtools::install(quick = TRUE)
# devtools::install()
devtools::load_all()
devtools::document(roclets = c('rd', 'collate', 'namespace'))


devtools::check()

devtools::build_manual()
devtools::build()


#usethis::use_vignette("survstan", "Introduction to the R package survstan")

devtools::build_vignettes()
devtools::build_readme()
devtools::build_site()

# # locally building site
# pkgdown::build_site()





# usethis::use_news_md()
# usethis::use_cran_comments()
devtools::test()
devtools::spell_check()
devtools::run_examples()
devtools::check()

devtools::release_checks()
devtools::missing_s3()
devtools::check_win_devel()
devtools::check_win_release()


cran_prep <- rhub::check_for_cran()
devtools::check_rhub(email = "fndemarqui@est.ufmg.br")
#devtools::check_rhub(email = "fndemarqui@est.ufmg.br")


# devtools::release(check = TRUE)
devtools::submit_cran()
