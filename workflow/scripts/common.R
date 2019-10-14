library("tidyverse")

get_prefix_col <- function(prefix, col_names) {

    covariate <- snakemake@params[["covariate"]]
    # add standard prefix to covariate
    col <- str_c(prefix, covariate, sep = "_")

    levels <- read_tsv(snakemake@input[["samples"]]) %>%
                select( !!covariate ) %>%
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
  
  
