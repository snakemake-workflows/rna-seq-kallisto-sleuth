rule sort_spia:
    input:
        spia_table="results/tables/pathways/{model}.pathways.tsv",
    output:
        spia_table_sorted="results/tables/pathways/{model}.pathways_sorted.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/yte/sort_spia/{model}.log",
    script:
        "../scripts/sort_spia.py"


rule render_datavzrd_config_spia:
    input:
        template=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        vega_circle=workflow.source_path(
            "../resources/custom_vega_plots/circle_diagram_de_genes.json"
        ),
        spia_table="results/tables/pathways/{model}.pathways_sorted.tsv",
    output:
        "results/datavzrd/spia/{model}.yaml",
    log:
        "logs/yte/render-datavzrd-config-spia/{model}.log",
    params:
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
    template_engine:
        "yte"


rule add_confidence_diffexp:
    input:
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        "results/tables/diffexp/{model}.genes-representative.diffexp_confidence.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/yte/add_confidence_diffexp/{model}.log",
    script:
        "../scripts/add_confidence_diffexp.py"


rule render_datavzrd_config_diffexp:
    input:
        template=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp_confidence.tsv",
        volcano_plots="results/plots/interactive/volcano/{model}.vl.json",
    output:
        "results/datavzrd/diffexp/{model}.yaml",
    params:
        samples=get_model_samples,
    log:
        "logs/yte/render-datavzrd-config-diffexp/{model}.log",
    template_engine:
        "yte"


rule copy_study_items_sig_terms:
    input:
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr{go_term_fdr}.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/yte/copy_study_items_sig_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/copy_study_items_sig_terms.py"


rule render_datavzrd_config_go_enrichment:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/go-enrichment-template.yaml"
        ),
        vega_circle=workflow.source_path(
            "../resources/custom_vega_plots/circle_diagram_genes.json"
        ),
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr{go_term_fdr}.tsv",
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


rule spia_datavzrd:
    input:
        config="results/datavzrd/spia/{model}.yaml",
        # files required for rendering the given configs
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
