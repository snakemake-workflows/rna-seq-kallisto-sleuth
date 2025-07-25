if is_3prime_experiment:
    kallisto_output = expand(
        "results/kallisto_3prime/{unit.sample}-{unit.unit}", unit=units.itertuples()
    )
else:
    kallisto_output = expand(
        "results/kallisto_cdna/{unit.sample}-{unit.unit}", unit=units.itertuples()
    )


rule compose_sample_sheet:
    input:
        config["samples"],
        config["units"],
        kallisto_output=kallisto_output,
    output:
        "results/sleuth/{model}.samples.tsv",
    log:
        "logs/{model}.compose-sample-sheet.log",
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
        samples="results/sleuth/{model}.samples.tsv",
        transcript_info="resources/transcripts_annotation.results.rds",
    output:
        sleuth_object="results/sleuth/{model,[^.]+}.rds",
        designmatrix="results/sleuth/{model}.designmatrix.rds",
    params:
        species=get_bioc_species_name(),
        model=get_model,
        exclude=config["diffexp"].get("exclude", None),
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/sleuth/{model}.init.log",
    group:
        "sleuth-init"
    threads: 6
    resources:
        # based on: https://github.com/pachterlab/sleuth/issues/139#issuecomment-331157007
        mem_mb=lambda wc, threads: threads * 8000,
    script:
        "../scripts/sleuth-init.R"


rule sleuth_diffexp:
    input:
        "results/sleuth/{model}.rds",
    output:
        mean_var_plot=report(
            "results/plots/mean-var/{model}.mean-variance-plot.pdf",
            caption="../report/plot-mean-var.rst",
            category="quality control",
            subcategory="per-model",
            labels={"model": "{model}", "plot": "mean-variance"},
        ),
        volcano_plots=report(
            "results/plots/volcano/{model}.volcano-plots.pdf",
            caption="../report/plot-volcano.rst",
            category="quality control",
            subcategory="per-model",
            labels={"model": "{model}", "plot": "volcano-plot"},
        ),
        ma_plots=report(
            "results/plots/ma/{model}.ma-plots.pdf",
            caption="../report/plot-ma.rst",
            category="quality control",
            subcategory="per-model",
            labels={"model": "{model}", "plot": "ma-plot"},
        ),
        qq_plots=report(
            "results/plots/qq/{model}.qq-plots.pdf",
            caption="../report/plot-qq.rst",
            category="quality control",
            subcategory="per-model",
            labels={"model": "{model}", "plot": "qq-plot"},
        ),
        transcripts_rds="results/sleuth/diffexp/{model}.transcripts.diffexp.rds",
        genes_aggregated_rds=(
            "results/sleuth/diffexp/{model}.genes-aggregated.diffexp.rds"
        ),
        genes_representative_rds=(
            "results/sleuth/diffexp/{model}.genes-representative.diffexp.rds"
        ),
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        sig_level_ma=config["diffexp"]["sig-level"]["ma-plot"],
        sig_level_qq=config["diffexp"]["sig-level"]["qq-plot"],
        representative_transcripts=config["resources"]["ref"][
            "representative_transcripts"
        ],
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
        results="results/tables/ihw/{model}.{level}.ihw-results.tsv",
        dispersion="results/plots/ihw/{level}/{model}.{level}.plot-dispersion.pdf",
        histograms="results/plots/ihw/{level}/{model}.{level}.plot-histograms.pdf",
        trends="results/plots/ihw/{level}/{model}.{level}.plot-trends.pdf",
        decision="results/plots/ihw/{level}/{model}.{level}.plot-decision.pdf",
        adj_pvals="results/plots/ihw/{level}/{model}.{level}.plot-adj-pvals.pdf",
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
            labels={"model": "{gene}-{transcript}-{model}"},
        ),
    conda:
        "../envs/sleuth.yaml"
    params:
        color_by=config["bootstrap_plots"]["color_by"],
        fdr=config["bootstrap_plots"]["FDR"],
        top_n=config["bootstrap_plots"]["top_n"],
        genes=config["diffexp"]["genes_of_interest"],
    log:
        "logs/plots/bootstrap/{model}/{model}.plot_bootstrap.log",
    script:
        "../scripts/plot-bootstrap.R"


rule prepare_pca:
    input:
        rds="results/sleuth/all.rds",
    output:
        # Write tsv instead of plot in order to create interactive plot with python since we did not find a good way to do it with R
        pca="results/plots/pca/{covariate}.pca.tsv",
        pc_var=report(
            "results/plots/pc-variance/{covariate}.pc-variance-plot.pdf",
            caption="../report/plot-pc-variance.rst",
            category="PCA",
            labels={"covariate": "{covariate}", "plot": "pc-variance-plot"},
        ),
        loadings=report(
            "results/plots/loadings/{covariate}.loadings-plot.pdf",
            caption="../report/plot-loadings.rst",
            category="PCA",
            labels={"covariate": "{covariate}", "plot": "loadings-plot"},
        ),
    conda:
        "../envs/sleuth.yaml"
    params:
        exclude_nas=config["pca"].get("exclude_nas", False),
    log:
        "logs/plots/pca/{covariate}.prepare_pca.log",
    script:
        "../scripts/prepare-pca.R"


rule plot_pca:
    input:
        pca="results/plots/pca/{covariate}.pca.tsv",
    output:
        pca=report(
            "results/plots/pca/{covariate}.pca.html",
            caption="../report/plot-pca.rst",
            category="PCA",
            labels={"covariate": "{covariate}", "plot": "pca"},
        ),
    conda:
        "../envs/pystats.yaml"
    log:
        "logs/plots/pca/{covariate}.plot_pca.log",
    params:
        color_by=lambda wildcards: wildcards.covariate,
    script:
        "../scripts/plot-pca.py"


rule plot_diffexp_pval_hist:
    input:
        diffexp_rds="results/sleuth/diffexp/{model}.{level}.diffexp.rds",
    output:
        report(
            "results/plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf",
            caption="../report/plot-pval-hist.rst",
            category="quality control",
            labels={
                "model": "{model}",
                "level": "{level}",
                "plot": "diffexp-pval-hist",
            },
        ),
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/plots/diffexp/{model}.{level}.diffexp-pval-hist.log",
    script:
        "../scripts/plot-diffexp-pval-hist.R"


rule logcount_matrix:
    input:
        "results/sleuth/{model}.rds",
    output:
        "results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/tables/logcount-matrix/{model}.logcount-matrix.log",
    script:
        "../scripts/sleuth-to-matrix.R"


rule tpm_matrix:
    input:
        "results/sleuth/{model}.rds",
    output:
        "results/tables/tpm-matrix/{model}.tpm-matrix.tsv",
    params:
        model=get_model,
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/tables/tpm-matrix/{model}.tpm-matrix.log",
    script:
        "../scripts/sleuth-to-tpm-matrix.R"


rule plot_diffexp_heatmap:
    input:
        logcountmatrix_file="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        predef_genelist=input_genelist,
    output:
        diffexp_heatmap=report(
            "results/plots/diffexp-heatmap/{model}.diffexp-heatmap.{mode}.pdf",
            caption="../report/plot-heatmap.rst",
            category="Heatmaps",
            labels={"model": "{model}-{mode}"},
        ),
    params:
        model=get_model,
    log:
        "logs/plots/diffexp-heatmap/{model}.diffexp-heatmap.{mode}.log",
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/plot_diffexp_heatmap.R"


rule plot_group_density:
    input:
        "results/sleuth/all.rds",
    output:
        report(
            "results/plots/group_density/{model}.group_density.pdf",
            caption="../report/plot-group-density.rst",
            category="quality control",
            labels={"model": "{model}-group_density"},
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
            category="quality control",
            labels={"model": "{model}-scatter-plot"},
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
            category="quality control",
            subcategory="per-sample",
            labels={"sample": "{sample}-{unit}", "plot": "fragment lengths"},
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
            category="quality control",
            labels={"model": "{model}-transcripts-plot-vars"},
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


rule vega_volcano_plot:
    input:
        tsv="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        spec=workflow.source_path("../../resources/vega_volcano_plot.json"),
    output:
        json="results/plots/interactive/volcano/{model}.vl.json",
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        primary_variable=lambda wc: config["diffexp"]["models"][wc.model][
            "primary_variable"
        ],
    log:
        "logs/vega-plots/volcano/{model}.log",
    conda:
        "../envs/vega.yaml"
    script:
        "../scripts/vega_plot_volcano.py"
