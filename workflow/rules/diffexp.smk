
kallisto_output = expand(
    "analysis/kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples())


rule compose_sample_sheet:
    input:
        kallisto_output,
        report(config["samples"], caption="../report/samples.rst")
    output:
        "analysis/sleuth/samples.tsv"
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
        return {"full": None}
    return config["diffexp"]["models"][wildcards.model]


rule sleuth_init:
    input:
        kallisto=kallisto_output,
        samples="analysis/sleuth/samples.tsv"
    output:
        "analysis/sleuth/{model,[^.]+}.rds"
    params:
        species=config["resources"]["ref"]["species"],
        model=lambda w: get_model(w)["full"],
        exclude=config["diffexp"].get("exclude", None)
    conda:
        "../envs/sleuth.yaml"
    group: "sleuth-init"
    script:
        "../scripts/sleuth-init.R"


checkpoint sleuth_diffexp:
    input:
        "analysis/sleuth/{model}.rds"
    output:
        transcripts_rds="analysis/sleuth/diffexp/{model}.transcripts.diffexp.rds",
        genes_aggregated_rds="analysis/sleuth/diffexp/{model}.genes-aggregated.diffexp.rds",
        genes_mostsigtrans_rds="analysis/sleuth/diffexp/{model}.genes-mostsigtrans.diffexp.rds",
        transcripts=report("analysis/tables/diffexp/{model}.transcripts.diffexp.tsv",
                            caption="../report/diffexp.rst",
                            category="Differential transcript expression"),
        genes_aggregated=report("analysis/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
                                caption="../report/diffexp-genes.rst",
                                category="Differential gene expression"),
        genes_mostsigtrans=report("analysis/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
                                    caption="../report/diffexp-mostsigtrans.rst",
                                    category="Differential gene expression")
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/sleuth-diffexp.R"


rule plot_bootstrap:
    input:
        "analysis/sleuth/{model}.rds"
    output:
        report("analysis/plots/bootstrap/{gene}/{gene}.{transcript}.{model}.bootstrap.pdf", caption="../report/plot-bootstrap.rst", category="Expression Plots")
    conda:
        "../envs/sleuth.yaml"
    params:
        color_by=config["bootstrap_plots"]["color_by"]
    script:
        "../scripts/plot-bootstrap.R"


rule plot_pca:
    input:
        "analysis/sleuth/all.rds"
    output:
        report("analysis/plots/pca/{covariate}.pca.pdf", caption="../report/plot-pca.rst", category="PCA")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-pca.R"


rule plot_diffexp_heatmap:
    input:
        so="analysis/sleuth/{model}.rds",
        diffexp="analysis/tables/diffexp/{model}.transcripts.diffexp.tsv"
    output:
        report("analysis/plots/diffexp-heatmap/{model}.diffexp-heatmap.pdf", caption="../report/heatmap.rst", category="Heatmaps")
    params:
        model=get_model
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-diffexp-heatmap.R"


rule plot_diffexp_pval_hist:
    input:
        diffexp_rds="analysis/sleuth/diffexp/{model}.{level}.diffexp.rds"
    output:
        report("analysis/plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf", caption="../report/pval-hist.rst", category="QC")
    params:
        model=get_model
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-diffexp-pval-hist.R"


rule tpm_matrix:
    input:
        "analysis/sleuth/{model}.rds"
    output:
        report("analysis/tables/tpm-matrix/{model}.tpm-matrix.tsv", caption="../report/tpm-matrix.rst", category="Expression Matrices")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/sleuth-to-matrix.R"


rule plot_fragment_length_dist:
    input:
        "analysis/sleuth/all.rds"
    output:
        report("analysis/plots/fld/{sample}-{unit}.fragment-length-dist.pdf", caption="../report/fld.rst", category="Fragment length distribution")
    conda:
        "../envs/sleuth.yaml"
    script:
        "../scripts/plot-fld.R"