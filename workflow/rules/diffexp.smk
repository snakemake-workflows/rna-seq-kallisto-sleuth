
kallisto_output = expand(
    "results/kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples())


rule compose_sample_sheet:
    input:
        kallisto_output,
        report(config["samples"], caption="../report/samples.rst")
    output:
        "results/sleuth/samples.tsv"
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
        samples="results/sleuth/samples.tsv"
    output:
        "results/sleuth/{model,[^.]+}.rds"
    params:
        species=config["resources"]["ref"]["species"],
        model=lambda w: get_model(w)["full"],
        exclude=config["diffexp"].get("exclude", None)
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/sleuth/{model}.init.log"
    group: "sleuth-init"
    threads: 6
    script:
        "../scripts/sleuth-init.R"


checkpoint sleuth_diffexp:
    input:
        "results/sleuth/{model}.rds"
    output:
        mean_var_plot=report("results/plots/mean-var/{model}.mean-variance-plot.pdf",
                            caption="../report/plot-mean-var.rst",
                            category="QC"),
        volcano_plots=report("results/plots/volcano/{model}.volcano-plots.pdf",
                            caption="../report/plot-volcano.rst",
                            category="QC"),
        ma_plots=report("results/plots/ma/{model}.ma-plots.pdf",
                            caption="../report/plot-ma.rst",
                            category="QC"),
        qq_plots=report("results/plots/qq/{model}.qq-plots.pdf",
                            caption="../report/plot-qq.rst",
                            category="QC"),
        transcripts_rds="results/sleuth/diffexp/{model}.transcripts.diffexp.rds",
        genes_aggregated_rds="results/sleuth/diffexp/{model}.genes-aggregated.diffexp.rds",
        genes_mostsigtrans_rds="results/sleuth/diffexp/{model}.genes-mostsigtrans.diffexp.rds",
        transcripts=report("results/tables/diffexp/{model}.transcripts.diffexp.tsv",
                            caption="../report/diffexp.rst",
                            category="Differential transcript expression"),
        genes_aggregated=report("results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
                                caption="../report/diffexp-genes.rst",
                                category="Differential gene expression"),
        genes_mostsigtrans=report("results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
                                    caption="../report/diffexp-mostsigtrans.rst",
                                    category="Differential gene expression")
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        sig_level_ma=config["diffexp"]["sig-level"]["ma-plot"],
        sig_level_qq=config["diffexp"]["sig-level"]["qq-plot"]
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/sleuth/{model}.diffexp.log"
    script:
        "../scripts/sleuth-diffexp.R"


rule plot_bootstrap:
    input:
        "results/sleuth/{model}.rds"
    output:
        report("results/plots/bootstrap/{gene}/{gene}.{transcript}.{model}.bootstrap.pdf", caption="../report/plot-bootstrap.rst", category="Expression Plots")
    conda:
        "../envs/sleuth.yaml"
    params:
        color_by=config["bootstrap_plots"]["color_by"]
    log:
        "logs/plots/bootstrap/{gene}.{transcript}.{model}.plot_bootstrap.log"
    script:
        "../scripts/plot-bootstrap.R"


rule plot_pca:
    input:
        "results/sleuth/all.rds"
    output:
        pca=report("results/plots/pca/{covariate}.pca.pdf", caption="../report/plot-pca.rst", category="PCA"),
        pc_var=report("results/plots/pc-variance/{covariate}.pc-variance-plot.pdf", caption="../report/pc-variance-plot.rst", category="PCA"),
        loadings=report("results/plots/loadings/{covariate}.loadings-plot.pdf", caption="../report/loadings-plot.rst", category="PCA")
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/pca/{covariate}.plot_pca.log"
    script:
        "../scripts/plot-pca.R"


rule plot_diffexp_heatmap:
    input:
        so="results/sleuth/{model}.rds",
        diffexp="results/tables/diffexp/{model}.transcripts.diffexp.tsv"
    output:
        report("results/plots/diffexp-heatmap/{model}.diffexp-heatmap.pdf", caption="../report/heatmap.rst", category="Heatmaps")
    params:
        model=get_model
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/diffexp-heatmap/{model}.diffexp-heatmap.log"
    script:
        "../scripts/plot-diffexp-heatmap.R"


rule plot_diffexp_pval_hist:
    input:
        diffexp_rds="results/sleuth/diffexp/{model}.{level}.diffexp.rds"
    output:
        report("results/plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf", caption="../report/pval-hist.rst", category="QC")
    params:
        model=get_model
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/diffexp/{model}.{level}.diffexp-pval-hist.log"
    script:
        "../scripts/plot-diffexp-pval-hist.R"


rule tpm_matrix:
    input:
        "results/sleuth/{model}.rds"
    output:
        report("results/tables/tpm-matrix/{model}.tpm-matrix.tsv", caption="../report/tpm-matrix.rst", category="Expression Matrices")
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/tables/tpm-matrix/{model}.tpm-matrix.log"
    script:
        "../scripts/sleuth-to-matrix.R"

rule plot_group_density:
    input:
        "results/sleuth/all.rds"
    output:
        report("results/plots/group_density/{model}.group_density.pdf", caption="../report/group-density.rst", category="QC")
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/group_density/{model}.group_density.log"
    script:
        "../scripts/plot-group-density.R"

rule plot_scatter:
     input:
         "results/sleuth/all.rds"
     output:
         report("results/plots/scatter/{model}.scatter.pdf", caption="../report/scatter.rst", category="QC")
     # params:
     #     covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
     conda:
         "../envs/sleuth.yaml"
     log:
         "logs/plots/scatter/{model}.scatter.log"
     script:
         "../scripts/plot-scatter.R"

rule plot_fragment_length_dist:
    input:
        "results/sleuth/all.rds"
    output:
        report("results/plots/fld/{sample}-{unit}.fragment-length-dist.pdf", caption="../report/fld.rst", category="Fragment length distribution")
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/fld/{sample}-{unit}.fragment-length-dist.log"
    script:
        "../scripts/plot-fld.R"

rule plot_vars:
    input:
        "results/sleuth/diffexp/{model}.transcripts.diffexp.rds"
    output:
        report("results/plots/variance/{model}.transcripts.plot_vars.pdf", caption="../report/plot-vars.rst", category="QC")
    params:
        model = get_model,
        sig_level = config["plot_vars"]["sig_level"]
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/variance/{model}.plot_vars.log"
    script:
        "../scripts/plot-variances.R"
