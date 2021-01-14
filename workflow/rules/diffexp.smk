kallisto_output = expand(
    "results/kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples()
)


rule compose_sample_sheet:
    input:
        report(config["samples"], caption="../report/samples.rst"),
        config["units"],
        kallisto_output=kallisto_output,
    output:
        "results/sleuth/samples.tsv",
    log:
        "logs/compose-sample-sheet.log",
    params:
        units=units,
        samples=samples,
    group:
        "sleuth-init"
    script:
        "../scripts/compose-sample-sheet.py"


rule sleuth_init:
    input:
        kallisto=kallisto_output,
        samples="results/sleuth/samples.tsv",
    output:
        "results/sleuth/{model,[^.]+}.rds",
    params:
        species=config["resources"]["ref"]["species"],
        model=lambda w: get_model(w)["full"],
        exclude=config["diffexp"].get("exclude", None),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/sleuth/{model}.init.log",
    group:
        "sleuth-init"
    threads: 6
    script:
        "../scripts/sleuth-init.R"


rule sleuth_diffexp:
    input:
        "results/sleuth/{model}.rds",
    output:
        mean_var_plot=report(
            "results/plots/mean-var/{model}.mean-variance-plot.pdf",
            caption="../report/plot-mean-var.rst",
            category="QC",
        ),
        volcano_plots=report(
            "results/plots/volcano/{model}.volcano-plots.pdf",
            caption="../report/plot-volcano.rst",
            category="QC",
        ),
        ma_plots=report(
            "results/plots/ma/{model}.ma-plots.pdf",
            caption="../report/plot-ma.rst",
            category="QC",
        ),
        qq_plots=report(
            "results/plots/qq/{model}.qq-plots.pdf",
            caption="../report/plot-qq.rst",
            category="QC",
        ),
        transcripts_rds="results/sleuth/diffexp/{model}.transcripts.diffexp.rds",
        genes_aggregated_rds=(
            "results/sleuth/diffexp/{model}.genes-aggregated.diffexp.rds"
        ),
        genes_mostsigtrans_rds=(
            "results/sleuth/diffexp/{model}.genes-mostsigtrans.diffexp.rds"
        ),
        transcripts=report(
            "results/tables/diffexp/{model}.transcripts.diffexp.tsv",
            caption="../report/diffexp.rst",
            category="Differential transcript expression",
        ),
        genes_aggregated=report(
            "results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
            caption="../report/diffexp-genes.rst",
            category="Differential gene expression",
        ),
        genes_mostsigtrans=report(
            "results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
            caption="../report/diffexp-mostsigtrans.rst",
            category="Differential gene expression",
        ),
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        sig_level_ma=config["diffexp"]["sig-level"]["ma-plot"],
        sig_level_qq=config["diffexp"]["sig-level"]["qq-plot"],
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/sleuth/{model}.diffexp.log",
    script:
        "../scripts/sleuth-diffexp.R"


rule ihw_fdr_control:
    input:
        "results/tables/diffexp/{model}.{level}.diffexp.tsv",
    output:
        results=report(
            "results/tables/ihw/{model}.{level}.ihw-results.tsv",
            caption="../report/ihw-results.rst",
            category="IHW",
        ),
        dispersion=report(
            "results/plots/ihw/{level}/{model}.{level}.plot-dispersion.pdf",
            caption="../report/plot-dispersion-ihw.rst",
            category="IHW",
        ),
        histograms=report(
            "results/plots/ihw/{level}/{model}.{level}.plot-histograms.pdf",
            caption="../report/plot-histograms-ihw.rst",
            category="IHW",
        ),
        trends=report(
            "results/plots/ihw/{level}/{model}.{level}.plot-trends.pdf",
            caption="../report/plot-trends-ihw.rst",
            category="IHW",
        ),
        decision=report(
            "results/plots/ihw/{level}/{model}.{level}.plot-decision.pdf",
            caption="../report/plot-decision-ihw.rst",
            category="IHW",
        ),
        adj_pvals=report(
            "results/plots/ihw/{level}/{model}.{level}.plot-adj-pvals.pdf",
            caption="../report/plot-adj-pvals-ihw.rst",
            category="IHW",
        ),
    conda:
        "../envs/ihw.yaml"
    log:
        "logs/tables/ihw/{model}.{level}.ihw.log",
    script:
        "../scripts/ihw-fdr-control.R"


rule plot_bootstrap:
    input:
        so="results/sleuth/{model}.rds",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
    output:
        report(
            directory("results/plots/bootstrap/{model}"),
            patterns=["{gene}.{transcript}.{model}.bootstrap.pdf"],
            caption="../report/plot-bootstrap.rst",
            category="Expression Plots",
        ),
    conda:
        "../envs/sleuth.yaml"
    params:
        color_by=config["bootstrap_plots"]["color_by"],
        fdr=config["bootstrap_plots"]["FDR"],
        top_n=config["bootstrap_plots"]["top_n"],
        genes=config["bootstrap_plots"]["genes_of_interest"],
    log:
        "logs/plots/bootstrap/{model}/{model}.plot_bootstrap.log",
    script:
        "../scripts/plot-bootstrap.R"


rule plot_pca:
    input:
        "results/sleuth/all.rds",
    output:
        pca=report(
            "results/plots/pca/{covariate}.pca.pdf",
            caption="../report/plot-pca.rst",
            category="PCA",
        ),
        pc_var=report(
            "results/plots/pc-variance/{covariate}.pc-variance-plot.pdf",
            caption="../report/plot-pc-variance.rst",
            category="PCA",
        ),
        loadings=report(
            "results/plots/loadings/{covariate}.loadings-plot.pdf",
            caption="../report/plot-loadings.rst",
            category="PCA",
        ),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/pca/{covariate}.plot_pca.log",
    script:
        "../scripts/plot-pca.R"


rule plot_diffexp_heatmap:
    input:
        so="results/sleuth/{model}.rds",
        diffexp="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
    output:
        report(
            "results/plots/diffexp-heatmap/{model}.diffexp-heatmap.pdf",
            caption="../report/plot-heatmap.rst",
            category="Heatmaps",
        ),
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/diffexp-heatmap/{model}.diffexp-heatmap.log",
    script:
        "../scripts/plot-diffexp-heatmap.R"


rule plot_diffexp_pval_hist:
    input:
        diffexp_rds="results/sleuth/diffexp/{model}.{level}.diffexp.rds",
    output:
        report(
            "results/plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf",
            caption="../report/plot-pval-hist.rst",
            category="QC",
        ),
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/diffexp/{model}.{level}.diffexp-pval-hist.log",
    script:
        "../scripts/plot-diffexp-pval-hist.R"


rule tpm_matrix:
    input:
        "results/sleuth/{model}.rds",
    output:
        report(
            "results/tables/tpm-matrix/{model}.tpm-matrix.tsv",
            caption="../report/tpm-matrix.rst",
            category="Expression Matrices",
        ),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/tables/tpm-matrix/{model}.tpm-matrix.log",
    script:
        "../scripts/sleuth-to-matrix.R"


rule plot_group_density:
    input:
        "results/sleuth/all.rds",
    output:
        report(
            "results/plots/group_density/{model}.group_density.pdf",
            caption="../report/plot-group-density.rst",
            category="QC",
        ),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/group_density/{model}.group_density.log",
    script:
        "../scripts/plot-group-density.R"


rule plot_scatter:
    input:
        "results/sleuth/all.rds",
    output:
        report(
            "results/plots/scatter/{model}.scatter.pdf",
            caption="../report/plot-scatter.rst",
            category="QC",
        ),
    # params:
    #     covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/scatter/{model}.scatter.log",
    script:
        "../scripts/plot-scatter.R"


rule plot_fragment_length_dist:
    input:
        "results/sleuth/all.rds",
    output:
        report(
            "results/plots/fld/{sample}-{unit}.fragment-length-dist.pdf",
            caption="../report/plot-fld.rst",
            category="Fragment length distribution",
        ),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/fld/{sample}-{unit}.fragment-length-dist.log",
    script:
        "../scripts/plot-fld.R"


rule plot_vars:
    input:
        "results/sleuth/diffexp/{model}.transcripts.diffexp.rds",
    output:
        report(
            "results/plots/variance/{model}.transcripts.plot_vars.pdf",
            caption="../report/plot-vars.rst",
            category="QC",
        ),
    params:
        model=get_model,
        sig_level=config["plot_vars"]["sig_level"],
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/variance/{model}.plot_vars.log",
    script:
        "../scripts/plot-variances.R"
