library("tidyverse")

load_bioconductor_package <- function(path_to_bioc_pkg_desc, pkg_name) {

    lib <- str_remove(path_to_bioc_pkg_desc, str_c(pkg_name, "/DESCRIPTION") )
    library("AnnotationDbi", lib.loc = lib)
    library(pkg_name, lib.loc = lib, character.only = TRUE)
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
