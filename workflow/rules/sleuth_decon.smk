## this rule is used to get a expression matrix on gene-level in which all 
## rows equal zero are removed and the normalized counts are log2(x + 3) 
## transformed

rule deconvolute_input_matrix:
    input:
        rds = "results/sleuth/null_model.rds",
    output:
        gene_counts = "results/tables/deconvolute/input_matrix.tsv",
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/tables/deconvolute/input_matrix.log",
    script:
        "../scripts/sleuth_decon.R"