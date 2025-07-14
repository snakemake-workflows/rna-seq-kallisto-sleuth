# General settings

To configure this workflow, modify the following files to reflect your dataset and differential expression analysis model:
* `config/samples.tsv`: samples sheet with covariates and conditions
* `config/units.tsv`: (sequencing) units sheet with raw data paths
* `config/config.yaml`: general workflow configuration and differential expression model setup


## samples sheet

For each biological sample, add a line to the sample sheet in `config/samples.tsv`.
The column `sample` is required and gives the sample name.
Additional columns can specify covariates (including batch effects) and conditions.
Please ensure that these column names can be interpreted as proper variable names.
Especially, they **should not** start with numeric chars (`1vs2`) or special (non-letter) symbols (`+compound_x`).
These columns can then be used in the `diffexp: models:` specification section in `config/config.yaml` (see below).

Missing values can be specified by empty columns or by writing `NA`.


## units sheet

For each sample, add one or more sequencing unit lines (runs, lanes or replicates) to the unit sheet in `config/units.tsv`.
Missing values can be specified by empty columns or by writing `NA`.

### input file options

For each unit, provide **one** out of the following options for input files:
* The path to two pairead-end FASTQ files in the columns `fq1`, `fq2`.
* The path to a single-end FASTQ file in the column `fq1`.
  For single-end data, you also need to specify `fragment_len_mean` and `fragment_len_sd`, which should usually be available from your sequencing provider.
* The path to a single-end BAM file in the column `bam_single`
* The path to a paired-end BAM file in the column `bam_paired`

### adapter trimming and read filtering

Finally, you can provide settings for the adapter trimming with `fastp` (see the [`fastp` documentation](https://github.com/OpenGene/fastp)):

In the column `fastp_adapters`, you can specify [known adapter sequences to be trimmed off by `fastp`](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters), including the command-line argument for the trimming.
For example, specify the following string in this column: `--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
If you don't know the adapters used, leave this empty (an empty string, containing no whitespace), and `fastp` will auto-detect the adapters that need to be trimmed.

In the column `fastp_extra`, you can specify [further `fastp` command-line settings](https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options).
If you leave this empty (an empty string, containing no whitespace), the workflow will set its default:
* [`--length_required 33`](https://github.com/OpenGene/fastp?tab=readme-ov-file#length-filter): This length filtering ensures that the resulting reads are all usable with the recommended k-mer size of `30` for `kallisto` quantification.
* [`--trim_poly_x --poly_x_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-X trimming removes polyA tails if they are 7 nucleotides or longer.
  It is run after adapter trimming.
* [`--trim_poly_g --poly_g_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-G trimming removes erroneous G basecalls at the tails of reads, if they are 7 nucleotides or longer.
  These Gs are artifacts in Illumina data from [machines with a one channel or two channel color chemistry](https://github.com/OpenGene/fastp/pull/508#issuecomment-3028078859).
  We currently set this by default, because [the auto-detection for the respective machines is lacking the latest machine types](https://github.com/OpenGene/fastp/pull/508).
  When the above-linked pull request is updated and merged, we can remove this and rely on the auto-detection.
If you want to specify additional command line options, we recommend always including those paramters in your units.tsv, as well.

#### Lexogen 3' QuantSeq adapter trimming

For this data, adapter trimming should automatically work as expected with the use of `fastp`.
The above-listed defaults are equivalent to an adaptation of the [Lexogen read preprocessing recommendations for 3' FWD QuantSeq data with `cutadapt`](https://faqs.lexogen.com/faq/what-sequences-should-be-trimmed).
The `fastp` equivalents, including minimal deviations from the recommendations, are motivated as follows:
* `-m`: In cutadapt, this is the short version of `--minimum-length`. The `fastp` equivalent is `--length_required`. In addition, we increase the minimum length to a default of 33 that makes more sense for kallisto quantification.
* `-O`: Here, `fastp` doesn't have an equivalent option, so we currently have to live with the suboptimal default of `4`. This is greater than the `min_overlap=3` used here; but smaller than the value of `7`, a threshold that we have found avoids removing randomly matching sequences when combined with the typical Illumina `max_error_rate=0.005`.
* `-a "polyA=A{20}"`: This can be replaced by with `fastp`'s dedicated option for `--trim_poly_x` tail removal ([which is run after adapter trimming](https://github.com/OpenGene/fastp?tab=readme-ov-file#global-trimming)).
* `-a "QUALITY=G{20}"`: This can be replaced by `fastp`'s dedicated option for the removal of artifactual trailing `G`s in Illumina data from [machines with a one channel or two channel color chemistry](https://github.com/OpenGene/fastp/pull/508#issuecomment-3028078859): `--trim_poly_g`.
  This is automatically activated for earlier Illumina machine models with this chemistry, but we recommend to activate it manually in the `fastp_extra` column of your `config/units.tsv` file for now, as [there are newer models that are not auto-detected, yet](https://github.com/OpenGene/fastp/pull/508).
  Also, we recommend to set `--poly_g_min_len 7`, to avoid trimming spurious matches of G-only stretches at the end of reads.
* `-n`: With the dedicated `fastp` options [getting applied in the right order](https://github.com/OpenGene/fastp?tab=readme-ov-file#global-trimming), this option is not needed any more.
* `-a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"`: We remove A{18}, as this is handled by `--trim_poly_x`.
  `fastp` uses the slightly higher `min_overlap` equivalent of `4`, which is currently hard-coded (and not exposed as a command-line argument).
  Because of this, we cannot set the `max_error_rate` to the Illumina error rate of about `0.005`.
* `-g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"`: This is not needed any more, as `fastp` searches the read sequence for adapter sequences from the start of the read (see [the `fastp` adapter search code](https://github.com/OpenGene/fastp/blob/723a4293a42f1ca05b93c37db6a157b4235c4dcc/src/adaptertrimmer.cpp#L92)).
* `--discard-trimmed`: We omit this, as adapter sequence removal from early on in the read will lead to short remaining read sequences being filtered by the `--length_required` argument.


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

For **transcript quantification**, `kallisto` is used.
For details regarding its command line arguments, see the [`kallisto` documentation](https://pachterlab.github.io/kallisto/manual).

#### Lexogen 3' QuantSeq data analysis

For Lexogen 3' QuantSeq data analysis, please set `experiment: 3-prime-rna-seq: activate: true` in the `config/config.yaml` file.
For more information information on Lexogen QuantSeq 3' sequencing, see: https://www.lexogen.com/quantseq-3mrna-sequencing/

### meta comparisons
Meta comparisons allow for comparing two full models against each other.
The axes represent the log2-fold changes (beta-scores) for the two models, with each point representing a gene. 
Points on the diagonal indicate no difference between the comparisons, while deviations from the diagonal suggest differences in gene expression between the treatments.
For more details see the comments in the `config.yaml`.
