def get_gmt_file(wildcards):
    return config["singscore"]["gene_sets"][wildcards.gene_set_group]["file"]

def get_color_col(wildcards):
    return config["singscore"]["gene_sets"][wildcards.gene_set_group]["plot_color_by"]

def get_shape_col(wildcards):
    return config["singscore"]["gene_sets"][wildcards.gene_set_group]["plot_shape_by"]

GENE_SET_GROUPS = list(config["singscore"]["gene_sets"].keys())

rule singscore:
    input:
        gene_counts = "results/tables/deconvolute/input_matrix.tsv",
        samples = config["samples"],
        gene_sets = get_gmt_file,
    output:
        "results/plots/singscore/{gene_set_group}.tiff"
    params:
        color_aes = get_color_col,
        shape_aes = get_shape_col,
    conda:
        "../envs/singscore.yaml"
    log:
        "logs/singscore_{gene_set_group}.log"
    script:
        "../scripts/singscore.R"