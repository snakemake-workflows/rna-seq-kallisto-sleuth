library("tidyverse")

load_bioconductor_package <- function(path_to_bioc_pkg, pkg_name) {

    lib <- str_remove(path_to_bioc_pkg, pkg_name)

    # ensure that dependencies of the pkg are also found at same location
    .libPaths( c( lib , .libPaths() ) )

    library(pkg_name, character.only = TRUE)

    print(str_c("loaded package ", pkg_name))

    # ensure that library() calls outside this function don't go looking in the
    # location needed here
    .libPaths( .libPaths()[-1] )
}

get_prefix_col <- function(prefix, col_names) {

    covariate <- snakemake@params[["covariate"]]
    # add standard prefix to covariate
    col <- str_c(prefix, covariate, sep = "_")

    levels <- read_tsv(snakemake@input[["samples"]]) %>%
                dplyr::select( !!covariate ) %>%
                distinct( ) %>%
                pull( !!covariate )

    # possible suffixes, cumulatively added on
    for (level in levels) {
        test_col <- col
        suffixes <- c("", level, ".0")
        found <- FALSE
        for(suffix in suffixes) {
            test_col <- str_c(test_col, suffix)
            # at the shortest hit possible, we break
            if(test_col %in% col_names) {
                return(test_col)
            }
        }
    }

    if(!found) {
        stop(str_c("Invalid covariate '", covariate, "', not found in diffexp table."))
    }

}
