from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path

##### load config and sample sheets #####


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False
)
samples.index.names = ["sample_id"]


def drop_unique_cols(df):
    singular_cols = df.nunique().loc[(df.nunique().values <= 1)].index
    return df.drop(singular_cols, axis=1)


samples = drop_unique_cols(samples)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(
    ["sample", "unit"], drop=False
)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


report: "../report/workflow.rst"


##### wildcard constraints #####


wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    model="|".join(list(config["diffexp"].get("models", [])) + ["all"]),


####### helpers ###########

is_3prime_experiment = (
    config.get("experiment", dict())
    .get("3-prime-rna-seq", dict())
    .get("activate", False)
)
three_prime_vendor = (
    config.get("experiment", dict()).get("3-prime-rna-seq", dict()).get("vendor")
)

if is_3prime_experiment:
    if three_prime_vendor != "lexogen":
        raise ValueError(
            f"Currently, only lexogene is supported. Please check the vendor "
            "in the config file and try again"
        )


def check_config():
    representative_transcripts_keywords = ["canonical", "mostsignificant"]
    representative_transcripts = config["resources"]["ref"][
        "representative_transcripts"
    ]
    if representative_transcripts not in representative_transcripts_keywords:
        if not os.path.exists(representative_transcripts):
            raise ValueError(
                f"Invalid value given for resources/ref/representative_transcripts in "
                "configuration. Must be 'canonical', 'mostsignificant' or valid path, "
                "but {representative_transcripts} does not exist or is not readable."
            )


check_config()


def is_activated(config_element):
    return config_element["activate"] in {"true", "True"}


def get_model(wildcards):
    if wildcards.model == "all":
        return {"full": None}
    return config["diffexp"]["models"][wildcards.model]


def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    fq2_not_present = pd.isnull(units.loc[(sample, unit), "fq2"])
    return fq2_not_present


def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[(wildcards.sample, wildcards.unit), "fq1"]
    else:
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        return [f"{u.fq1}", f"{u.fq2}"]


def get_all_fastqs(wildcards):
    for item in units[["sample", "unit"]].itertuples():
        if is_single_end(item.sample, item.unit):
            yield f"results/trimmed/{item.sample}-{item.unit}.fastq.gz"
        else:
            yield f"results/trimmed/{item.sample}-{item.unit}.1.fastq.gz"
            yield f"results/trimmed/{item.sample}-{item.unit}.2.fastq.gz"


def get_model_samples(wildcards):
    samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    units = pd.read_csv(config["units"], sep="\t", dtype=str, comment="#")
    sample_file = units.merge(samples, on="sample")
    sample_file["sample_name"] = sample_file.apply(
        lambda row: "{}-{}".format(row["sample"], row["unit"]), axis=1
    )
    gps = config["diffexp"]["models"][wildcards.model]["primary_variable"]
    sample_groups = sample_file.loc[sample_file[gps].notnull(), ["sample_name"]]
    samples = sample_groups["sample_name"].values
    return samples


def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards,
        )
    # single end sample
    return expand("results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)


def get_bioc_species_name():
    first_letter = config["resources"]["ref"]["species"][0]
    subspecies = config["resources"]["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_bioc_species_pkg():
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = get_bioc_species_name()[0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)


def render_enrichment_env():
    species_pkg = f"bioconductor-{get_bioc_species_pkg()}"
    with open(workflow.source_path("../envs/enrichment.yaml")) as f:
        env = yaml.load(f, Loader=yaml.SafeLoader)
    env["dependencies"].append(species_pkg)
    env_path = Path("resources/envs/enrichment.yaml")
    env_path.parent.mkdir(parents=True, exist_ok=True)
    with open(env_path, "w") as f:
        yaml.dump(env, f)
    return env_path.absolute()


bioc_species_pkg = get_bioc_species_pkg()
enrichment_env = render_enrichment_env()


def kallisto_quant_input(wildcards):
    if is_3prime_experiment:
        return "results/canonical_reads/{sample}-{unit}.fastq"
    elif not is_single_end(wildcards.sample, wildcards.unit):
        return expand(
            "results/trimmed/{{sample}}-{{unit}}.{group}.fastq.gz", group=[1, 2]
        )
    else:
        return expand("results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)


def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if len(input.fastq) == 1 or is_3prime_experiment:
        extra += " --single --single-overhang --pseudobam"
        extra += (
            " --fragment-length {unit.fragment_len_mean} " "--sd {unit.fragment_len_sd}"
        ).format(unit=units.loc[(wildcards.sample, wildcards.unit)])
    else:
        extra += " --fusion"
    return extra


def input_genelist(predef_genelist):
    if config["diffexp"]["genes_of_interest"]["activate"] == True:
        predef_genelist = config["diffexp"]["genes_of_interest"]["genelist"]
    else:
        predef_genelist = []

    return predef_genelist


def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """

    wanted_input = []

    # request goatools if 'activated' in config.yaml
    if config["enrichment"]["goatools"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
                    "results/plots/go_terms/{model}.go_term_enrichment_{go_ns}.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.pdf",
                    "results/datavzrd-reports/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}",
                ],
                model=config["diffexp"]["models"],
                go_ns=["BP", "CC", "MF"],
                gene_fdr=str(config["enrichment"]["goatools"]["fdr_genes"]).replace(
                    ".", "-"
                ),
                go_term_fdr=str(
                    config["enrichment"]["goatools"]["fdr_go_terms"]
                ).replace(".", "-"),
            )
        )

    # request fgsea if 'activated' in config.yaml
    if config["enrichment"]["fgsea"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/fgsea/{model}.all-gene-sets.tsv",
                    "results/tables/fgsea/{model}.sig-gene-sets.tsv",
                    "results/plots/fgsea/{model}.table-plot.pdf",
                    "results/plots/fgsea/{model}",
                ],
                model=config["diffexp"]["models"],
            )
        )
    # request spia if 'activated' in config.yaml
    if config["enrichment"]["spia"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/pathways/{model}.pathways.tsv",
                    "results/datavzrd-reports/spia-{model}/",
                ],
                model=config["diffexp"]["models"],
            )
        )

    # workflow output that is always wanted
    # general sleuth output
    wanted_input.extend(
        expand(
            [
                "results/plots/mean-var/{model}.mean-variance-plot.pdf",
                "results/plots/volcano/{model}.volcano-plots.pdf",
                "results/plots/ma/{model}.ma-plots.pdf",
                "results/plots/qq/{model}.qq-plots.pdf",
                "results/tables/diffexp/{model}.transcripts.diffexp.tsv",
                "results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
                "results/sleuth/{model}.samples.tsv",
                "results/datavzrd-reports/diffexp-{model}",
                "results/plots/diffexp-heatmap/{model}.diffexp-heatmap.{mode}.pdf",
            ],
            model=config["diffexp"]["models"],
            mode=["topn"],
        )
    )
    if config["diffexp"]["genes_of_interest"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/plots/diffexp-heatmap/{model}.diffexp-heatmap.{mode}.pdf",
                ],
                model=config["diffexp"]["models"],
                mode=["predefined"],
            )
        )

    # ihw false discovery rate control
    wanted_input.extend(
        expand(
            [
                "results/tables/ihw/{model}.{level}.ihw-results.tsv",
                "results/plots/ihw/{level}/{model}.{level}.plot-dispersion.pdf",
                "results/plots/ihw/{level}/{model}.{level}.plot-histograms.pdf",
                "results/plots/ihw/{level}/{model}.{level}.plot-trends.pdf",
                "results/plots/ihw/{level}/{model}.{level}.plot-decision.pdf",
                "results/plots/ihw/{level}/{model}.{level}.plot-adj-pvals.pdf",
            ],
            model=config["diffexp"]["models"],
            level=["transcripts", "genes-aggregated", "genes-representative"],
        )
    )

    # sleuth p-value histogram plots
    wanted_input.extend(
        expand(
            "results/plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf",
            model=config["diffexp"]["models"],
            level=["transcripts", "genes-aggregated", "genes-representative"],
        )
    )

    # technical variance vs. observed variance
    # wanted_input.extend(
    #        expand("results/plots/variance/{model}.transcripts.plot_vars.pdf", model=config["diffexp"]["models"]),
    #    )

    # PCA plots of kallisto results, each coloured for a different covariate
    wanted_input.extend(
        expand(
            [
                "results/plots/pc-variance/{covariate}.pc-variance-plot.pdf",
                "results/plots/loadings/{covariate}.loadings-plot.pdf",
                "results/plots/pca/{covariate}.pca.pdf",
            ],
            covariate=samples.columns[samples.columns != "sample"],
        )
    )

    # group-density plot
    wanted_input.extend(
        expand(
            ["results/plots/group_density/{model}.group_density.pdf"],
            model=config["diffexp"]["models"],
        )
    )

    # scatter plots
    if config["scatter"]["activate"]:
        wanted_input.extend(
            expand(
                ["results/plots/scatter/{model}.scatter.pdf"],
                model=config["diffexp"]["models"],
            )
        )

    # sleuth bootstrap plots
    wanted_input.extend(
        expand("results/plots/bootstrap/{model}", model=config["diffexp"]["models"])
    )

    # fragment length distribution plots
    wanted_input.extend(
        expand(
            "results/plots/fld/{unit.sample}-{unit.unit}.fragment-length-dist.pdf",
            unit=units[["sample", "unit"]].itertuples(),
        )
    )

    if is_activated(config["diffsplice"]):
        # diffsplice analysis
        wanted_input.extend(
            expand(
                "results/plots/diffsplice/{model}/{cons}",
                model=config["diffexp"]["models"],
                cons=["with_consequences", "without_consequences"],
            )
        )

    if is_3prime_experiment:
        wanted_input.extend(
            expand(
                "results/plots/QC/3prime-QC-plot.{ind_transcripts}.html",
                ind_transcripts=config["experiment"]["3-prime-rna-seq"]["plot-qc"],
            )
        )

    if (
        is_3prime_experiment
        and config["experiment"]["3-prime-rna-seq"]["plot-qc"] != "all"
    ):
        wanted_input.extend(
            expand(
                "results/plots/QC/3prime-ind-QC-plot.{ind_transcripts}.html",
                ind_transcripts=config["experiment"]["3-prime-rna-seq"]["plot-qc"],
            )
        )
    return wanted_input
