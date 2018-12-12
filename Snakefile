include: "rules/common.smk"
include: "rules/quant.smk"
include: "rules/diffexp.smk"

rule all:
    input:
        expand("tables/diffexp/{model}.diffexp.tsv", model=config["diffexp"]["models"]),
