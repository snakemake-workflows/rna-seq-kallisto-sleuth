# General settings

To configure this workflow, modify the following files to reflect your dataset and differential expression analysis model:
* `config/samples.tsv`: samples sheet with covariates and conditions
* `config/units.tsv`: (sequencing) units sheet with raw data paths
* `config/config.yaml`: general workflow configuration and differential expression model setup

# samples sheet

For each biological sample, add a line to the sample sheet in `config/samples.tsv`.
The column `sample` is required and gives the sample name.
Additional columns can specify covariates (including batch effects) and conditions.
These columns can then be used in the `diffexp: models:` specification section in `config/config.yaml` (see below)

Missing values can be specified by empty columns or by writing `NA`.

# units sheet

For each sample, add one or more sequencing unit lines (runs, lanes or replicates) to the unit sheet in `config/units.tsv`.
For each unit, provide either of the following:
* The path to two pairead-end FASTQ files in the columns `fq1`, `fq2`.
* The path to a single-end FASTQ file in the column `fq1`.
  For single-end data, you also need to specify `fragment_len_mean` and `fragment_len_sd`, which should usually be available from your sequencing provider.

Missing values can be specified by empty columns or by writing `NA`.

# config.yaml

This file contains the general workflow configuration and the setup for the differential expression analysis performed by sleuth.
Configurable options should be explained in the comments above the respective entry or right here in this `config/README.md` section.
If something is unclear, don't hesitate to [file an issue in the `rna-seq-kallisto-sleuth` GitHub repository](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/new/choose).

## differential expression model setup

The core functionality of this workflow is provided by the software [`sleuth`](https://pachterlab.github.io/sleuth/about).
You can use it to test for differential expression of genes or transcripts between two or more subgroups of samples.

### main sleuth model

The main idea of sleuth's internal model, is to test a `full:` model (containing (a) variable(s) of interest AND batch effects) against a `reduced:` model (containing ONLY the batch effects).
So these are the most important entries to set up under any model that you specify via `diffexp: models:`.
If you don't know any batch effects, the `reduced:` model will have to be `~1`.
Otherwise it will be the tilde followed by an addition of the names of any columns that contain batch effects, for example: `reduced: ~batch_effect_1 + batch_effect_2`.
The full model than additionally includes variables of interest, so fore example: `full: ~variable_of_interest + batch_effect_1 + batch_effect_2`.

### sleuth effect sizes

Effect size estimates are calculated as so-called beta-values by `sleuth`.
For binary comparisons (your variable of interest has two factor levels), they resemble a log2 fold change.
To know which variable of interest to use for the effect size calculation, you need to provide its column name as the `primary_variable:`.
And for sleuth to know what level of that variable of interest to use as the base level, specify the respective entry as the `base_level:`.
