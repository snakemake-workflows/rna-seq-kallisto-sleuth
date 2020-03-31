**Volcano plots** computed with sleuth for wald test using the model ``{{ snakemake.params.model["full"] }}``.
The plots display beta values (regression coefficient) on the x-axis vs. log(significance) on the y-axis and has a significance level of {{ snakemake.params.sig_level_volcano }}.
Significant genes are coloured red.
