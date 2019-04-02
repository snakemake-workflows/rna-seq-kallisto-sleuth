from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


####### helpers ###########

def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    return units.loc[
        (wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return expand("trimmed/{sample}-{unit}.fastq.gz", **wildcards)

def get_bootstrap_plots(wildcards):
    """Dynamically determine which transcripts to plot based on
       checkpoint output."""
    transcripts = dict()
    genes = set()
    for model in config["diffexp"]["models"]:
        # Obtain results from the sleuth_diffexp checkpoint.
        # This happens dynamically after the checkpoint is completed, and
        # is skipped automatically before completion.
        results = pd.read_csv(
            checkpoints.sleuth_diffexp.get(model=model).output[0], sep="\t")
        # group transcripts by gene
        genes.update(results[results.qval <= config["diffexp"]["FDR"]].ext_gene)
        for g in genes:
            trx = set()
            trx.update(results[results.ext_gene == g][results.qval <= config["diffexp"]["FDR"]].target_id)
            transcripts[g] = trx
    # Require the respective output from the plot_bootstrap rule.
    return ["plots/bootstrap/{gene}/{transcript}.bootstrap.svg".format(gene=g, transcript=t)
            for t in transcripts[g] for g in genes]
