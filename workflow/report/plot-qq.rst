**Q-Q plots** computed with sleuth using the model ``{{ snakemake.params.model["full"] }}``.
Plots show theoretical quantile expected from a standard normal distribution on x-axis vs. observed quantiles from wald test and likelihood ratio test.
For likelihood ratio test are used aggregated p-values.
For wald test the plots are created for ``{{ snakemake.wildcards.covariate }}`` and have an additional blue colored approximation line which passes through the first and third quartile.
Significant genes up to significance level of ``{{ snakemake.params.sig_level_qq }}`` are coloured red.
