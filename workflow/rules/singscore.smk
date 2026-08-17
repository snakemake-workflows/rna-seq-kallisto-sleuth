rule singscore:
    input:
        gene_counts="results/tables/tpm-matrix/{model}.tpm-matrix.sorted.tsv",
        samples=lookup(within=config, dpath="samples"),
        gene_sets=lookup(
            within=config, dpath="enrichment/singscore/gene_sets/{gene_set_group}/file"
        ),
        transcript_ref="resources/transcripts_annotation.results.tsv",
    output:
        "results/plots/singscore/{model}.{gene_set_group}.tiff",
    log:
        "logs/singscore_{model}.{gene_set_group}.log",
    conda:
        "../envs/singscore.yaml"
    params:
        up_set=lookup(
            within=config, 
            dpath="enrichment/singscore/gene_sets/{gene_set_group}/upregulated_gene_set",
            ),
        down_set=lookup(
            within=config, 
            dpath="enrichment/singscore/gene_sets/{gene_set_group}/downregulated_gene_set",
            ),
        color_aes=lookup(
            within=config,
            dpath="enrichment/singscore/gene_sets/{gene_set_group}/plot_color_by",
        ),
        shape_aes=lookup(
            within=config,
            dpath="enrichment/singscore/gene_sets/{gene_set_group}/plot_shape_by",
        ),
    script:
        "../scripts/singscore.R"
