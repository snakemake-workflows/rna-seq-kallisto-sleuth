**MA plots** computed with sleuth for wald test using the model ``{{ snakemake.params.model["full"] }}``.
The plots display, for each transcript, the mean of abundances across samples on the x-axis and fold change on the y-axis.
Significant genes up to significance level of {{ snakemake.params.sig_level_ma }} are coloured red.
