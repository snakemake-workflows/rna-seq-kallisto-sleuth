**Gene set enrichment results table** as determined by fgsea (with nperm=``{{ snakemake.params.nperm }}`` random gene sets drawn from the background of all measured genes for each gene set tested, to generate an empirical null distribution to test against). The most significant transcript of each gene was determined by sleuth using the model ``{{ snakemake.params.model["full"] }}``.

