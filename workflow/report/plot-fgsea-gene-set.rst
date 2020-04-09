**Gene set enrichment plot** of gene set ``{{ snakemake.wildcards.gene_set }}``, found to be significantly enriched for differentially expressed genes by fgsea.
The most significant transcript of each gene was determined by sleuth using the model ``{{ snakemake.params.model["full"] }}``.

