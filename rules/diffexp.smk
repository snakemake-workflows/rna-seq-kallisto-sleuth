
kallisto_output = expand(
    "kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples())


rule compose_sample_sheet:
    input:
        kallisto_output
    output:
        "sleuth/samples.tsv"
    group: "sleuth-init"
    run:
        samples_ = units[["sample", "unit"]].merge(samples, on="sample")
        samples_["sample"] = samples_.apply(
            lambda row: "{}-{}".format(row["sample"], row["unit"]), axis=1)
        samples_["path"] = kallisto_output
        del samples_["unit"]
        samples_.to_csv(output[0], sep="\t")


rule sleuth_init:
    input:
        kallisto=kallisto_output,
        samples="sleuth/samples.tsv"
    output:
        "sleuth/all.rds"
    params:
        species=config["ref"]["species"]
    conda:
        "../envs/sleuth.yaml"
    group: "sleuth-init"
    script:
        "../scripts/sleuth-init.R"


checkpoint sleuth_diffexp:
    input:
        "sleuth/all.rds"
    output:
        report("tables/diffexp/{model}.diffexp.tsv", caption="../report/diffexp.rst"),
        report("tables/diffexp/{model}.aggregated.diffexp.tsv", caption="../report/diffexp.rst")
    params:
        model=lambda wildcards: config["diffexp"]["models"][wildcards.model]
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/sleuth-diffexp.R"


rule plot_bootstrap:
    input:
        "sleuth/all.rds"
    output:
        report("plots/bootstrap/{gene}/{transcript}.bootstrap.svg", caption="../report/plot-bootstrap.rst")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-bootstrap.R"


rule plot_pca:
    input:
        "sleuth/all.rds"
    output:
        report("plots/pca/{covariate}.pca.svg", caption="../report/plot-pca.rst")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-pca.R"
