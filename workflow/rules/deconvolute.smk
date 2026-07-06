## this rule is used to make deconvolution calculation based on the  
## normalized log2(x + 3) transformed counts produced by the rule sleuth_decon

species = config["resources"]["ref"]["species"]

rule deconvolution:
    input:
        gene_counts = "results/tables/deconvolute/input_matrix.tsv",
        ref = f"resources/xcell2/{species}/{{ref_set}}.rda",
    output:
        f"results/tables/deconvolute/plots/heatmap/{{ref_set}}.tiff"
    conda:
        "../envs/deconvolute.yaml"
    log:
        "logs/tables/deconvolute/input_matrix.log", # no idea yet
    script:
        "../scripts/deconvolute.R"