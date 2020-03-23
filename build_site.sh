R -e "devtools::document()"
R -e "pkgdown::build_site(document = FALSE, examples = FALSE, preview = FALSE)"
