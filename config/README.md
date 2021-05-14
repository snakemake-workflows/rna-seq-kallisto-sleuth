# General settings

To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample and unit sheet

* Add samples to `config/samples.tsv`. For each sample, the column `sample` and some `condition` column need to be specified.
  If you want to name your `condition` differently, please also adjust the `primary_variable` of the model specification in the `config/config.yaml` accordingly.
  Similarly, you can include columns for `batch_effect`s and have to specify them accordingly in the `full` and `reduced` versions of your model in the `config/config.yaml`.
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`.
  For each unit, provide either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). 
  If only one FASTQ file is provided (single-end data), you also need to specify `fragment_len_mean` and `fragment_len_sd`, which should usually be available from your sequencing provider.

Missing values can be specified by empty columns or by writing `NA`.


