include: "rules/common.smk"
include: "rules/quant.smk"

rule all:
    input:
        expand("kallisto/{u.sample}-{u.unit}", u=units.itertuples())
