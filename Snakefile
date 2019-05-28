include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/quant.smk"
include: "rules/diffexp.smk"

rule all:
    input:
        expand("tables/diffexp/{model}.diffexp.tsv", model=config["diffexp"]["models"]),
        expand("plots/heatmap/{model}.heatmap.pdf", model=config["diffexp"]["models"]),
        expand("plots/pca/{covariate}.pca.pdf", covariate=samples.columns[samples.columns != "sample"]),
        [get_bootstrap_plots(model) for model in config["diffexp"]["models"]]
