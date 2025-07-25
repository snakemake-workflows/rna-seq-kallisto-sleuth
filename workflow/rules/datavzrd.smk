# Postprocessing GO Enrichment Data
rule postprocess_go_enrichment:
    input:
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
    conda:
        "../envs/polars.yaml"
    log:
        "logs/yte/postprocess_go_enrichment/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/postprocess_go_enrichment.py"


# Postprocessing Differential Expression Data
# Does not work for level = genes-aggregated since it does not contain beta values.
rule postprocess_diffexp:
    input:
        genes_representative="results/tables/diffexp/{model}.{level}.diffexp.tsv",
    output:
        "results/tables/diffexp/{model}.{level}.diffexp_postprocessed.tsv",
    conda:
        "../envs/pandas.yaml"
    params:
        model=get_model,
    log:
        "logs/yte/postprocess_diffexp/{model}/{level}.log",
    script:
        "../scripts/postprocess_diffexp.py"


rule postprocess_tpm_matrix:
    input:
        tpm="results/tables/tpm-matrix/{model}.tpm-matrix.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp_postprocessed.tsv",
    output:
        "results/tables/tpm-matrix/{model}.tpm-matrix.sorted.tsv",
    conda:
        "../envs/pandas.yaml"
    params:
        model=get_model,
    log:
        "logs/tables/tpm-matrix/{model}.tpm-matrix.sort.log",
    script:
        "../scripts/postprocess_tpm.py"


rule plot_enrichment_scatter:
    input:
        "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
    output:
        "results/plots/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}_scatter.json",
    conda:
        "../envs/pystats.yaml"
    params:
        identifier="term",
        effect_x="effect",
        effect_y="p_fdr_bh",
    log:
        "logs/plot_enrichment_scatter-{model}/plot_enrichment_scatter-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/plot_enrichment_pathway_scatter.py"


rule plot_pathway_scatter:
    input:
        "results/tables/pathways/{model}.{database}.pathways.tsv",
    output:
        "results/plots/pathways/{model}.{database}.pathways.tsv_scatter.json",
    conda:
        "../envs/pystats.yaml"
    log:
        "logs/plot_pathway_scatter/{model}.{database}.pathways.tsv_scatter.log",
    params:
        identifier="Name",
        effect_x="total perturbation accumulation",
        effect_y="Combined FDR",
    script:
        "../scripts/plot_enrichment_pathway_scatter.py"


# Generating SPIA Datavzrd Report
rule spia_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        # files required for rendering the given configs
        vega_circle=workflow.source_path(
            "../resources/custom_vega_plots/circle_diagram_genes.json"
        ),
        spia_table="results/tables/pathways/{model}.{database}.pathways.tsv",
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
        scatter="results/plots/pathways/{model}.{database}.pathways.tsv_scatter.json",
    output:
        report(
            directory("results/datavzrd-reports/spia-{model}_{database}"),
            htmlindex="index.html",
            caption="../report/spia.rst",
            category="Pathway enrichment",
            patterns=["index.html"],
            labels={
                "model": "{model}",
                "database": "{database}",
            },
        ),
    log:
        "logs/datavzrd-report/spia-{model}/spia-{model}_{database}.log",
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
    wrapper:
        "v5.5.0/utils/datavzrd"


# Generating Differential Expression Datavzrd Report
rule diffexp_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        # optional files required for rendering the given config
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        tpm_matrix="results/tables/tpm-matrix/{model}.tpm-matrix.sorted.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp_postprocessed.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp_postprocessed.tsv",
        volcano_plots="results/plots/interactive/volcano/{model}.vl.json",
    output:
        report(
            directory("results/datavzrd-reports/diffexp-{model}"),
            htmlindex="index.html",
            caption="../report/diffexp.rst",
            category="Differential expression analysis",
            patterns=["index.html"],
            labels={"model": "{model}"},
        ),
    log:
        "logs/datavzrd-report/diffexp.{model}/diffexp.{model}.log",
    params:
        extra="",
        model=get_model,
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        samples=get_model_samples,
        primary_variable=lambda wildcards: config["diffexp"]["models"][
            wildcards.model
        ]["primary_variable"],
    wrapper:
        "v5.5.0/utils/datavzrd"


# Generating GO Enrichment Datavzrd Report
rule go_enrichment_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/go-enrichment-template.yaml"),
        vega_bars=workflow.source_path(
            "../resources/custom_vega_plots/horizontal_bars.json"
        ),
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
        scatter="results/plots/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}_scatter.json",
    output:
        report(
            directory(
                "results/datavzrd-reports/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}"
            ),
            htmlindex="index.html",
            caption="../report/go-enrichment-sig_terms.rst",
            category="GO term enrichment",
            subcategory="{model}",
            patterns=["index.html"],
            labels={
                "model": "{model}",
                "gene_fdr": "{gene_fdr}",
                "go_term_fdr": "{go_term_fdr}",
            },
        ),
    log:
        "logs/datavzrd-report/go_enrichment-{model}/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        samples=get_model_samples,
    wrapper:
        "v5.5.0/utils/datavzrd"


# Generating Meta Comparison Datavzrd Reports
rule meta_compare_datavzrd:
    input:
        config=lambda wildcards: workflow.source_path(
            f"../resources/datavzrd/meta_comparison-{wildcards.method}-template.yaml"
        ),
        table="results/tables/{method}/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/{method}/{meta_comp}.json",
    output:
        report(
            directory("results/datavzrd-reports/{method}_meta_comparison_{meta_comp}"),
            htmlindex="index.html",
            caption="../report/meta_compare.rst",
            category="Comparisons",
            subcategory="{meta_comp}",
            patterns=["index.html"],
            labels=lambda wildcards: get_meta_compare_labels(
                method=f"{wildcards.method.capitalize()}: "
            )(wildcards),
        ),
    log:
        "logs/datavzrd-report/meta_comp_{method}.{meta_comp}.log",
    wrapper:
        "v5.5.0/utils/datavzrd"


# Generating Input Datavzrd Reports
rule inputs_datavzrd:
    input:
        config=lambda wc: workflow.source_path(
            f"../resources/datavzrd/{wc.input}-template.yaml"
        ),
        table=lambda wc: config[wc.input],
    output:
        report(
            directory("results/datavzrd-reports/inputs/{input}"),
            htmlindex="index.html",
            category="Inputs",
            patterns=["index.html"],
            labels={
                "input": "{input}",
            },
        ),
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
    log:
        "logs/datavzrd-report/{input}_datavzrd.log",
    wrapper:
        "v5.5.0/utils/datavzrd"
