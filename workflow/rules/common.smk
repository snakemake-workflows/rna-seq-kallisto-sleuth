from snakemake.utils import validate
from itertools import product
import pandas as pd


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")
report: "../report/workflow.rst"

##### wildcard constraints #####

wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"])


####### helpers ###########

def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    fq2_present = pd.isnull(units.loc[(sample, unit), "fq2"])
    if isinstance(fq2_present, pd.core.series.Series):
        # if this is the case, get_fastqs cannot work properly
        raise ValueError(
            f"Multiple fq2 entries found for sample-unit combination {sample}-{unit}.\n"
            "This is most likely due to a faulty units.tsv file, e.g. "
            "a unit name is used twice for the same sample.\n"
            "Try checking your units.tsv for duplicates."
        )
    return fq2_present

### check for each sample, that...
for s, r in units.groupby("sample"):
    if not ( # all units are single end
            all(map(lambda x: is_single_end(x[0], x[1]), product(s, r.unit.values))) or
            # all units are paired end
            all(map(lambda x: not is_single_end(x[0], x[1]), product(s, r.unit.values)))
            ):
        raise ValueError("kallisto requires units within a sample to either\n"
                            "all be paired end, or all be single end.\n"
                            f"Sample {s} has a mix, please fix.")

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def get_trimmed(wildcards):
    files=[]
    sample=wildcards.sample
    us=units.loc[sample, "unit"].tolist()
    for unit in us:
        if not is_single_end(sample, unit):
            # paired-end sample
            files.extend(
                expand( [ "results/trimmed/{sample}-{unit}.{group}.fastq.gz" ],
                        sample=sample, unit=unit,  group=[1, 2])
                        )
        else:
            files.extend([ f"results/trimmed/{sample}-{unit}.fastq.gz" ])
    return files

def get_bioc_species_pkg(wildcards):
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = config["resources"]["ref"]["species"][0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)

def get_bioc_pkg_path(wildcards):
    return "resources/bioconductor/lib/R/library/{pkg}".format(pkg=get_bioc_species_pkg(wildcards))

def is_activated(config_element):
    return config_element['activate'] in {"true","True"}

def get_bootstrap_plots(model, gene_list=None):
    """Dynamically determine which transcripts to plot based on
       checkpoint output."""
    transcripts = dict()
    genes = set()
    # Obtain results from the sleuth_diffexp checkpoint.
    # This happens dynamically after the checkpoint is completed, and
    # is skipped automatically before completion.
    results = pd.read_csv(
        checkpoints.sleuth_diffexp.get(model=model).output.transcripts, sep="\t").dropna()
    # group transcripts by gene
    if gene_list is None:
        results = results[results.qval <= config["bootstrap_plots"]["FDR"]][:config["bootstrap_plots"]["top_n"]]
        genes = set(results.ext_gene)
    else:
        genes = set(gene_list)
        genes.add("Custom")
    for g in genes:
        if not pd.isnull(g):
            valid = results.ext_gene == g
            trx = set(results[valid].target_id)
            transcripts[g] = trx
    # Require the respective output from the plot_bootstrap rule.
    return ["results/plots/bootstrap/{gene}/{gene}.{transcript}.{model}.bootstrap.pdf".format(gene=g, transcript=t, model=model)
            for g, ts in transcripts.items()
            for t in ts]

def get_fgsea_plots(model):
    plots = set()
    table = pd.read_csv(
                checkpoints.fgsea.get( model=model ).output.significant,
                sep="\t").dropna()
    gs = set(table['pathway'])
    for gene_set in gs:
        plots.add(
            "results/plots/fgsea/{model}.{gene_set}.gene-set-plot.pdf".format(
                model=model,
                gene_set=gene_set
                )
            )
    return plots
