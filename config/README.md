# General settings

To configure this workflow, modify the following files to reflect your dataset and differential expression analysis model:
* `config/samples.tsv`: samples sheet with covariates and conditions
* `config/units.tsv`: (sequencing) units sheet with raw data paths
* `config/config.yaml`: general workflow configuration and differential expression model setup

## samples sheet

For each biological sample, add a line to the sample sheet in `config/samples.tsv`.
The column `sample` is required and gives the sample name.
Additional columns can specify covariates (including batch effects) and conditions.
These columns can then be used in the `diffexp: models:` specification section in `config/config.yaml` (see below)

Missing values can be specified by empty columns or by writing `NA`.

## units sheet

For each sample, add one or more sequencing unit lines (runs, lanes or replicates) to the unit sheet in `config/units.tsv`.
For each unit, provide either of the following:
* The path to two pairead-end FASTQ files in the columns `fq1`, `fq2`.
* The path to a single-end FASTQ file in the column `fq1`.
  For single-end data, you also need to specify `fragment_len_mean` and `fragment_len_sd`, which should usually be available from your sequencing provider.
* The path to a single-end BAM file in the column `bam_single`
* The path to a paired-end bam BAM file in the column `bam_paired`

Missing values can be specified by empty columns or by writing `NA`.

## config.yaml

This file contains the general workflow configuration and the setup for the differential expression analysis performed by sleuth.
Configurable options should be explained in the comments above the respective entry or right here in this `config/README.md` section.
If something is unclear, don't hesitate to [file an issue in the `rna-seq-kallisto-sleuth` GitHub repository](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/new/choose).

### differential expression model setup

The core functionality of this workflow is provided by the software [`sleuth`](https://pachterlab.github.io/sleuth/about).
You can use it to test for differential expression of genes or transcripts between two or more subgroups of samples.

#### main sleuth model

The main idea of sleuth's internal model, is to test a `full:` model (containing (a) variable(s) of interest AND batch effects) against a `reduced:` model (containing ONLY the batch effects).
So these are the most important entries to set up under any model that you specify via `diffexp: models:`.
If you don't know any batch effects, the `reduced:` model will have to be `~1`.
Otherwise it will be the tilde followed by an addition of the names of any columns that contain batch effects, for example: `reduced: ~batch_effect_1 + batch_effect_2`.
The full model than additionally includes variables of interest, so fore example: `full: ~variable_of_interest + batch_effect_1 + batch_effect_2`.

#### sleuth effect sizes

Effect size estimates are calculated as so-called beta-values by `sleuth`.
For binary comparisons (your variable of interest has two factor levels), they resemble a log2 fold change.
To know which variable of interest to use for the effect size calculation, you need to provide its column name as the `primary_variable:`.
And for sleuth to know what level of that variable of interest to use as the base level, specify the respective entry as the `base_level:`.

### preprocessing `params`

For **adapter trimming**, `cutadapt` is used, with some defaults for standard Illumina data given in the `config.yaml`.
For more details see the comments in the `config.yaml` or the [`cutadapt` documentation](https://cutadapt.readthedocs.io).
For parameter suggestions for Lexogen 3' QuantSeq data, see the section below.

For **transcript quantification**, `kallisto` is used.
For details regarding its command line arguments, see the [`kallisto` documentation](https://pachterlab.github.io/kallisto/manual).

#### Lexogen 3' QuantSeq data analysis

For Lexogen 3' QuantSeq data analysis, please set `experiment: 3-prime-rna-seq: activate: true` in the `config/config.yaml` file.
For more information information on Lexogen QuantSeq 3' sequencing, see: https://www.lexogen.com/quantseq-3mrna-sequencing/
In addition, for Lexogen 3' FWD QuantSeq data, we recommend setting the `params: cutadapt-se:` with:
```
    adapters: "-a 'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=7;max_error_rate=0.005'"
    extra: "--minimum-length 33 --nextseq-trim=20 --poly-a"
```
This is an adaptation of the [Lexogen read preprocessing recommendations for 3' FWD QuantSeq data](https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-1).
Changes to the recommendations are motivated as follows:
* `-m`: We switched to the easier to read `--minimum-length` and apply this minimum length globally. In addition, we increase the minimum length to a default of 33 that makes more sense for kallisto quantification.
* `-O`: Instead of this option, minimum overlap is specified per expression.
* `-a "polyA=A{20}"`: We replace this with `cutudapt`s dedicated option for `--poly-a` tail removal (which is run after adapter trimming).
* `-a "QUALITY=G{20}"`: We replace this with `cutudapt`s dedicated option for the removal artifactual trailing `G`s in NextSeq and NovaSeq data: `--nextseq-trim=20`.
* `-n`: With the dedicated `cutadapt` options getting applied in the right order, this option is not needed any more.
* `-a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"`: We remove A{18}, as this is handled by `--poly-a`. We increase `min_overlap` to 7 and set the `max_error_rate` to the Illumina error rate of about 0.005, both to avoid spurious adapter matches being removed.
* `-g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"`: This is not needed any more, as `-a` option will lead to complete removal of read sequence if adapter is found at the start of the read, see: https://cutadapt.readthedocs.io/en/stable/guide.html#rightmost
* `--discard-trimmed`: We omit this, as the `-a` with the adapter sequence will lead to complete read sequence removal if adapter is found at start, and the `--minimum-length` will then discard such empty reads.
