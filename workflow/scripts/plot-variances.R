log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")

sl <- snakemake@params[['sig_level']]

so <- sleuth_load(snakemake@input[[1]])

# colors from colorblind-safe Brewer palette "Dark2":
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
cb_safe_green <- "#1B9E77"
cb_safe_red <- "#D95F02"

p <- plot_vars(so,
        test = NULL,
        test_type = "lrt",
        which_model = "full",
        sig_level = sl,
        point_alpha = 0.2,
        sig_color = cb_safe_red,
        xy_line = TRUE,
        xy_line_color = cb_safe_red,
        highlight = NULL,
        highlight_color = cb_safe_green
    )
ggsave(filename = snakemake@output[[1]], width = 7, height = 7)
