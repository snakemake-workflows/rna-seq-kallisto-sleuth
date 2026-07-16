rule singscore:
    input:
        gene_counts = "results/tables/tpm-matrix/{model}.tpm-matrix.sorted.tsv",
        samples = lookup(within = config, dpath="samples"),
        gene_sets = lookup(within = config, dpath="enrichment/singscore/gene_sets/{gene_set_group}/file"),
    output:
        "results/plots/singscore/{gene_set_group}.tiff",
    params:
        color_aes = lookup(within = config, dpath="enrichment/singscore/gene_sets/{gene_set_group}/plot_color_by"),
        shape_aes = lookup(within = config, dpath="enrichment/singscore/gene_sets/{gene_set_group}/plot_shape_by"),
    conda:
        "../envs/singscore.yaml"
    log:
        "logs/singscore_{gene_set_group}.log"
    script:
        "../scripts/singscore.R"