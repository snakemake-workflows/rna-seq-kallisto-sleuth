library("tidyverse")

get_beta_col <- function(covariate, col_names) {

    # add standard prefix to covariate
    beta_col <- str_c("b", covariate, sep = "_")

    # possible suffixes, cumulatively added on
    suffixes <- c("", "1", ".0")
    found <- FALSE
    for(suffix in suffixes) {
        beta_col <- str_c(beta_col, suffix)
        # at the shortest hit possible, we break
        if(beta_col %in% col_names) {
            found <- TRUE
            break
        }
    }

    if(!found) {
        stop(str_c("Invalid covariate ", covariate, ", not found in diffexp table."))
    } else {
        return(beta_col)
    }

}
  
  
