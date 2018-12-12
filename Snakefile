include: "rules/common.smk"
include: "rules/quant.smk"
include: "rules/diffexp.smk"

rule all:
    input:
        "sleuth/all.rds"
