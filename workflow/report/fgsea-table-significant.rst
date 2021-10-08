**Significantly enriched gene sets** as determined by fgsea with a Benjamini-Hochberg controlled false discovery rate control below ``{{ snakemake.params.gene_set_fdr }}``. fgsea determined an empirical null distribution for each tested gene set from nperm=``{{ snakemake.params.nperm }}`` random gene sets drawn from the background of all genes measured.
Using the sleuth model ``{{ snakemake.params.model["full"] }}``.

