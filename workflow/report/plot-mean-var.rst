**Plot of mean-variance relationship of transcripts** computed with sleuth using the model ``{{ snakemake.params.model["full"] }}``.
The plot displays the mean expression of transcripts pooled across all samples on x-axis and the biological variance (after removing technical variance) on y-axis.
Transcripts that are used in the shrinkage estimation are colored blue.
The fitted curve is colored red and represents the smooth line.
