# Postprocessing GO Enrichment Data
rule postprocess_go_enrichment:
    input:
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/yte/postprocess_go_enrichment/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/postprocess_go_enrichment.py"


# Postprocessing Differential Expression Data
rule postprocess_diffexp:
    input:
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        "results/tables/diffexp/{model}.genes-representative.diffexp_postprocessed.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/yte/postprocess_diffexp/{model}.log",
    script:
        "../scripts/postprocess_diffexp.py"


# Rendering Datavzrd Config for Differential Expression
rule render_datavzrd_config_diffexp:
    input:
        template=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp_postprocessed.tsv",
        volcano_plots="results/plots/interactive/volcano/{model}.vl.json",
    output:
        "results/datavzrd/diffexp/{model}.yaml",
    params:
        samples=get_model_samples,
    log:
        "logs/yte/render-datavzrd-config-diffexp/{model}.log",
    template_engine:
        "yte"


# Rendering Datavzrd Config for SPIA Pathway Enrichment Analysis
rule render_datavzrd_config_spia:
    input:
        template=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        vega_circle=workflow.source_path(
            "../resources/custom_vega_plots/circle_diagram_de_genes.json"
        ),
        spia_table="results/tables/pathways/{model}.pathways.tsv",
    output:
        "results/datavzrd/spia/{model}.yaml",
    log:
        "logs/yte/render-datavzrd-config-spia/{model}.log",
    params:
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
    template_engine:
        "yte"


# Rendering Datavzrd Config for GO Enrichment
rule render_datavzrd_config_go_enrichment:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/go-enrichment-template.yaml"
        ),
        vega_bars=workflow.source_path(
            "../resources/custom_vega_plots/horizontal_bars.json"
        ),
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        "results/datavzrd/go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.yaml",
    params:
        samples=get_model_samples,
    log:
        "logs/yte/render-datavzrd-config-go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    template_engine:
        "yte"


# Generating SPIA Datavzrd Report
rule spia_datavzrd:
    input:
        config="results/datavzrd/spia/{model}.yaml",
        spia_table="results/tables/pathways/{model}.pathways.tsv",
    output:
        report(
            directory("results/datavzrd-reports/spia-{model}"),
            htmlindex="index.html",
            caption="../report/spia.rst",
            category="Pathway enrichment",
            patterns=["index.html"],
            labels={"model": "{model}"},
        ),
    log:
        "logs/datavzrd-report/spia-{model}/spia-{model}.log",
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"


# Generating Differential Expression Datavzrd Report
rule diffexp_datavzrd:
    input:
        config="results/datavzrd/diffexp/{model}.yaml",
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
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
    params:
        model=get_model,
    log:
        "logs/datavzrd-report/diffexp.{model}/diffexp.{model}.log",
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"


# Generating GO Enrichment Datavzrd Report
rule go_enrichment_datavzrd:
    input:
        config="results/datavzrd/go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.yaml",
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        sig_go="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        report(
            directory(
                "results/datavzrd-reports/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}"
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
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"


# Meta Comparison Configuration
rule render_datavzrd_config_meta_comparison:
    input:
        template=lambda wildcards: workflow.source_path(
            f"../resources/datavzrd/meta_comparison-{wildcards.method}-template.yaml"
        ),
        table="results/tables/{method}/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/{method}/{meta_comp}.json",
    output:
        "results/datavzrd/{method}/meta_comparison_{meta_comp}.yaml",
    log:
        "logs/yte/render-datavzrd-config-meta_comparison/{method}/{meta_comp}.log",
    template_engine:
        "yte"


# Generating Pathway Meta Comparison Datavzrd Report
rule meta_compare_pathways_datavzrd:
    input:
        config="results/datavzrd/pathways/meta_comparison_{meta_comp}.yaml",
        table="results/tables/pathways/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/pathways/{meta_comp}.json",
    output:
        report(
            directory("results/datavzrd-reports/pathways_meta_comparison_{meta_comp}"),
            htmlindex="index.html",
            caption="../report/meta_compare.rst",
            category="Comparisons",
            subcategory="{meta_comp}",
            patterns=["index.html"],
            labels=get_meta_compare_labels(method="Pathway: "),
        ),
    log:
        "logs/datavzrd-report/meta_comp_pathways.{meta_comp}.log",
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"


# Generating Differential Expression Meta Comparison Datavzrd Report
rule meta_compare_diffexp_datavzrd:
    input:
        config="results/datavzrd/diffexp/meta_comparison_{meta_comp}.yaml",
        table="results/tables/diffexp/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/diffexp/{meta_comp}.json",
    output:
        report(
            directory("results/datavzrd-reports/diffexp_meta_comparison_{meta_comp}"),
            htmlindex="index.html",
            caption="../report/meta_compare.rst",
            category="Comparisons",
            subcategory="{meta_comp}",
            patterns=["index.html"],
            labels=get_meta_compare_labels(method="Diffexp: "),
        ),
    log:
        "logs/datavzrd-report/meta_comp_diffexp.{meta_comp}.log",
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"


# Generating GO Terms Meta Comparison Datavzrd Report
rule meta_compare_enrichment_datavzrd:
    input:
        config="results/datavzrd/go_terms/meta_comparison_{meta_comp}.yaml",
        table="results/tables/go_terms/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/go_terms/{meta_comp}.json",
    output:
        report(
            directory("results/datavzrd-reports/go_terms_meta_comparison_{meta_comp}"),
            htmlindex="index.html",
            caption="../report/meta_compare.rst",
            category="Comparisons",
            subcategory="{meta_comp}",
            patterns=["index.html"],
            labels=get_meta_compare_labels("Enrichment: "),
        ),
    log:
        "logs/datavzrd-report/meta_comp_go_terms.{meta_comp}.log",
    wrapper:
        "v3.10.2-3-gbeb9d22/utils/datavzrd"
