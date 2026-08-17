if is_long_read_sequencing:
    kallisto_output = expand(
        "results/kallisto_long_cdna/{unit.sample}-{unit.unit}/quant-tcc",
        unit=units.itertuples(),
    )
else:
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
    group:
        "sleuth-init"
    params:
        units=units,
        samples=samples,
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
    log:
        "logs/sleuth/{model}.init.log",
    group:
        "sleuth-init"
    conda:
        "../envs/sleuth.yaml"
    threads: 6
    resources:
        # based on: https://github.com/pachterlab/sleuth/issues/139#issuecomment-331157007
        mem_mb=lambda wc, threads: threads * 8000,
    params:
        species=get_bioc_species_name(),
        model=get_model,
        exclude=config["diffexp"].get("exclude", None),
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
    log:
        "logs/sleuth/{model}.diffexp.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        sig_level_ma=config["diffexp"]["sig-level"]["ma-plot"],
        sig_level_qq=config["diffexp"]["sig-level"]["qq-plot"],
        representative_transcripts=config["resources"]["ref"][
            "representative_transcripts"
        ],
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
    log:
        "logs/tables/ihw/{model}.{level}.ihw.log",
    conda:
        "../envs/ihw.yaml"
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
    log:
        "logs/plots/bootstrap/{model}/{model}.plot_bootstrap.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        color_by=config["bootstrap_plots"]["color_by"],
        fdr=config["bootstrap_plots"]["FDR"],
        top_n=config["bootstrap_plots"]["top_n"],
        genes_of_interest=lookup(within=config, dpath="diffexp/genes_of_interest"),
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
    log:
        "logs/plots/pca/{covariate}.prepare_pca.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        exclude_nas=config["pca"].get("exclude_nas", False),
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
    log:
        "logs/plots/pca/{covariate}.plot_pca.log",
    conda:
        "../envs/pystats.yaml"
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
    log:
        "logs/plots/diffexp/{model}.{level}.diffexp-pval-hist.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        model=get_model,
    script:
        "../scripts/plot-diffexp-pval-hist.R"


rule logcount_matrix:
    input:
        "results/sleuth/{model}.rds",
    output:
        "results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
    log:
        "logs/tables/logcount-matrix/{model}.logcount-matrix.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        model=get_model,
    script:
        "../scripts/sleuth-to-matrix.R"


rule tpm_matrix:
    input:
        "results/sleuth/{model}.rds",
    output:
        "results/tables/tpm-matrix/{model}.tpm-matrix.tsv",
    log:
        "logs/tables/tpm-matrix/{model}.tpm-matrix.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        model=get_model,
    script:
        "../scripts/sleuth-to-tpm-matrix.R"


rule plot_diffexp_heatmap:
    input:
        logcountmatrix_file="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        # default provides a dummy file path for the "topn" case (a file path
        # that is always valid, but never loaded)
        predef_gene_list=lookup(
            within=config,
            dpath="diffexp/genes_of_interest/gene_lists/{gene_list}",
            default=lookup(within=config, dpath="samples"),
        ),
    output:
        diffexp_heatmap=report(
            "results/plots/diffexp-heatmap/{model}.diffexp-heatmap.{gene_list}.pdf",
            caption="../report/plot-heatmap.rst",
            category="Heatmaps",
            labels={
                "model": "{model}",
                "gene list": "{gene_list}",
            },
        ),
    log:
        "logs/plots/diffexp-heatmap/{model}.diffexp-heatmap.{gene_list}.log",
    conda:
        "../envs/heatmap.yaml"
    params:
        model=get_model,
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
    log:
        "logs/plots/group_density/{model}.group_density.log",
    conda:
        "../envs/sleuth.yaml"
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
    log:
        "logs/plots/scatter/{model}.scatter.log",
    # params:
    #     covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/sleuth.yaml"
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
    log:
        "logs/plots/fld/{sample}-{unit}.fragment-length-dist.log",
    conda:
        "../envs/sleuth.yaml"
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
    log:
        "logs/plots/variance/{model}.plot_vars.log",
    conda:
        "../envs/sleuth.yaml"
    params:
        model=get_model,
        sig_level=config["plot_vars"]["sig_level"],
    script:
        "../scripts/plot-variances.R"


rule vega_volcano_plot:
    input:
        tsv="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        spec=workflow.source_path("../../resources/vega_volcano_plot.json"),
    output:
        json="results/plots/interactive/volcano/{model}.vl.json",
    log:
        "logs/vega-plots/volcano/{model}.log",
    conda:
        "../envs/vega.yaml"
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        primary_variable=lambda wc: config["diffexp"]["models"][wc.model][
            "primary_variable"
        ],
    script:
        "../scripts/vega_plot_volcano.py"
