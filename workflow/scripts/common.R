suppressPackageStartupMessages({
    library("tidyverse")
})

get_prefix_col <- function(covariate, prefix, col_names) {

    # add standard prefix to covariate
    col <- str_c(prefix, covariate, sep = "_")

    # possible suffixes, cumulatively added on
    suffixes <- c("", "1", ".0")
    found <- FALSE
    for(suffix in suffixes) {
        col <- str_c(col, suffix)
        # at the shortest hit possible, we break
        if(col %in% col_names) {
            found <- TRUE
            break
        }
    }

    if(!found) {
        stop(str_c("Invalid covariate ", covariate, ", not found in diffexp table."))
    } else {
        return(col)
    }

}
  
  
