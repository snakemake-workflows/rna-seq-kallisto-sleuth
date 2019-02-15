include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/quant.smk"
include: "rules/diffexp.smk"

rule all:
    input:
        expand("tables/diffexp/{model}.diffexp.tsv", model=config["diffexp"]["models"]),
        expand("plots/pca/{covariate}.pca.svg", covariate=samples.columns[samples.columns != "sample"]),
        get_bootstrap_plots
