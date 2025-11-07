**Gene set enrichment results table** as determined by fgsea (with nperm=``{{ snakemake.params.nperm }}`` random gene sets drawn from the background of all measured genes for each gene set tested, to generate an empirical null distribution to test against), using the sleuth model ``{{ snakemake.params.model["full"] }}``.

