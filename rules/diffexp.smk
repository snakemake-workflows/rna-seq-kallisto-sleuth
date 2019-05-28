
kallisto_output = expand(
    "kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples())


rule compose_sample_sheet:
    input:
        kallisto_output,
        report(config["samples"], caption="../report/samples.rst")
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


def get_model(wildcards):
    if wildcards.model == "all":
        return None
    return config["diffexp"]["models"][wildcards.model]["full"]


rule sleuth_init:
    input:
        kallisto=kallisto_output,
        samples="sleuth/samples.tsv"
    output:
        "sleuth/{model}.rds"
    params:
        species=config["ref"]["species"],
        model=get_model
    conda:
        "../envs/sleuth.yaml"
    group: "sleuth-init"
    script:
        "../scripts/sleuth-init.R"


checkpoint sleuth_diffexp:
    input:
        "sleuth/{model}.rds"
    output:
        transcripts=report("tables/diffexp/{model}.diffexp.tsv", caption="../report/diffexp.rst"),
        genes=report("tables/diffexp/{model}.aggregated.diffexp.tsv", caption="../report/diffexp-genes.rst"),
        heatmap=report("plots/heatmap/{model}.heatmap.pdf", caption="../report/heatmap.rst")
    params:
        model=get_model,
        reduced_model=lambda wildcards: config["diffexp"]["models"][wildcards.model]["reduced"]
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/sleuth-diffexp.R"


rule plot_bootstrap:
    input:
        "sleuth/{model}.rds"
    output:
        report("plots/bootstrap/{gene}.{transcript}.{model}.bootstrap.pdf", caption="../report/plot-bootstrap.rst")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-bootstrap.R"


rule plot_pca:
    input:
        "sleuth/all.rds"
    output:
        report("plots/pca/{covariate}.pca.pdf", caption="../report/plot-pca.rst")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-pca.R"
