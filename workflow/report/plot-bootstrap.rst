Distribution of estimated expression of {{ snakemake.wildcards.transcript }} over all bootstrap runs of kallisto.
Values are given in transcripts per million (TPM) derived from normalized counts.
Bootstrap plots are generated for all ``genes_of_interest`` in the ``config.yaml`` file and for the {{ snakemake.params.top_n }} transcripts with the lowest qvalue among the transcripts with an FDR below {{ snakemake.params.fdr }} for the model ``{{ snakemake.wildcards.model }}``.
