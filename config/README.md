# General settings

To configure this workflow, modify the following files to reflect your dataset and differential expression analysis model:
* `config/samples.tsv`: samples sheet with covariates and conditions
* `config/units.tsv`: (sequencing) units sheet with raw data paths
* `config/config.yaml`: general workflow configuration and differential expression model setup

For the `samples.tsv` and `units.tsv`, we explain the expected colummns right here, in this `README.md` file.
For the `config.yaml` file, all entries are explained in detail in its comments.


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

For each unit, provide exactly **one** out of the following options for input files:
* The path to two paired-end FASTQ files in the columns `fq1`, `fq2`.
* The path to a single-end FASTQ file in the column `fq1`.
  For single-end data, you also need to specify `fragment_len_mean` and `fragment_len_sd`, which should usually be available from your sequencing provider.
* The path to a single-end BAM file in the column `bam_single`
* The path to a paired-end BAM file in the column `bam_paired`

### adapter trimming and read filtering

Finally, you can provide settings for the adapter trimming with `fastp` (see the [`fastp` documentation](https://github.com/OpenGene/fastp)) via the `units.tsv` columns `fastp_adapters` and `fastp_extra`.
However, if you leave those two columns empty (no whitespace!), `fastp` will auto-detect adapters and the workflow will set sensible defaults for trimming of RNA-seq data.
If you use this automatic inference, make sure to double-check the `Detected read[12] adapter:` entries in the resulting `fastp` HTML report.
This is part of the final `snakemake` report of the workflow, or can be found in the sample-specific folders under `results/trimmed/`, once a sample has been processed this far.
If the auto-detection didn't work at all (empty `Detected read[12] adapter:` entries), or the `Occurrences` in the `Adapters` section are lower than you would expect, please ensure that you find out which adapters were used and configure the adapter trimming manually:

In the column `fastp_adapters`, you can specify [known adapter sequences to be trimmed off by `fastp`](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters), including the command-line argument for the trimming.
For example, specify the following string in this column: `--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
If you don't know the adapters used, leave this empty (an empty string, containing no whitespace), and `fastp` will auto-detect the adapters that need to be trimmed.
If you want to make the auto-detection explicit for paired-end samples, you can also specify `--detect_adapter_for_pe`.

In the column `fastp_extra`, you can specify [further `fastp` command-line settings](https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options).
If you leave this empty (an empty string, containing no whitespace), the workflow will set its default:
* [`--length_required 33`](https://github.com/OpenGene/fastp?tab=readme-ov-file#length-filter): This length filtering ensures that the resulting reads are all usable with the recommended k-mer size of `30` for `kallisto` quantification.
* [`--trim_poly_x --poly_x_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-X trimming removes polyA tails if they are 7 nucleotides or longer.
  It is run after adapter trimming.
* [`--trim_poly_g --poly_g_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-G trimming removes erroneous G basecalls at the tails of reads, if they are 7 nucleotides or longer.
  These Gs are artifacts in Illumina data from [machines with a one channel or two channel color chemistry](https://github.com/OpenGene/fastp/pull/508#issuecomment-3028078859).
  We currently set this by default, because [the auto-detection for the respective machines is lacking the latest machine types](https://github.com/OpenGene/fastp/pull/508).
  When the above-linked pull request is updated and merged, we can remove this and rely on the auto-detection.
If you want to specify additional command line options, we recommend always including those parameters in your units.tsv, as well.
Here's the full concatenation for copy-pasting:

```bash
--length_required 33 --trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7
```

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
* `--discard-trimmed`: We omit this, as adapter sequence removal early in the read will leave short remaining read sequences that are subsequently filtered by `--length_required`.


## config.yaml

This file contains the general workflow configuration and the setup for the differential expression analysis performed by sleuth.
Configurable options should be explained in the comments above the respective entry, so the easiest way to set it up for your workflow is to carefully read through the `config/config.yaml` file and adjust it to your needs.
If something is unclear, don't hesitate to [file an issue in the `rna-seq-kallisto-sleuth` GitHub repository](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/new/choose).
