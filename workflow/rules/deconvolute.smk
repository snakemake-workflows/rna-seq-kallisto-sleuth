## this rule is used to get a expression matrix on gene-level in which all 
## rows equal zero are removed and the normalized counts are log2(x + 3) 
## transformed

# rule deconvolute_input_matrix:
#     input:
#         rds = "results/sleuth/null_model.rds",
#     output:
#         gene_counts = "results/tables/deconvolute/input_matrix.tsv",
#     conda:
#         "../envs/sleuth.yaml"
#     log:
#         "logs/tables/deconvolute/input_matrix.log",
#     script:
#         "../scripts/sleuth_decon.R"
        
## this rule is used to make deconvolution calculation based on the  
## normalized log2(x + 3) transformed counts produced by the rule sleuth_decon

rule deconvolution:
    input:
        gene_counts="results/tables/tpm-matrix/{model}.tpm-matrix.tsv",
        ref="resources/celltype_references/{celltype_reference}.xCell2Ref.rda",
    output:
        heatmap="results/plots/deconvolute/{celltype_reference}.{model}.heatmap.tiff",
        decon_scores="results/tables/deconvolute/{celltype_reference}.{model}.decon_scores.tsv",
    conda:
        "../envs/deconvolute.yaml"
    log:
        "logs/tables/deconvolute/{celltype_reference}.{model}.log"
    script:
        "../scripts/deconvolute.R"