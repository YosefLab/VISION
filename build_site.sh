R -e "devtools::document()"
R -e "pkgdown::build_site(devel = FALSE, examples = FALSE, preview = FALSE)"
