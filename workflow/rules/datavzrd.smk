rule render_datavzrd_config_spia:
    input:
        template=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        spia_table="results/tables/pathways/{model}.pathways.tsv",
        spia_table_activated="results/tables/pathways/{model}.activated-pathways.tsv",
        spia_table_inhibited="results/tables/pathways/{model}.inhibited-pathways.tsv",
    output:
        "results/datavzrd/spia/{model}.yaml",
    log:
        "logs/yte/render-datavzrd-config-spia/{model}.log",
    template_engine:
        "yte"


rule render_datavzrd_config_diffexp:
    input:
        template=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        volcano_plots="results/plots/interactive/volcano/{model}.vl.json",
    output:
        "results/datavzrd/diffexp/{model}.yaml",
    params:
        samples=get_model_samples,
    log:
        "logs/yte/render-datavzrd-config-diffexp/{model}.log",
    template_engine:
        "yte"


rule render_datavzrd_config_go_enrichment:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/go-enrichment-template.yaml"
        ),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
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
        spia_table_activated="results/tables/pathways/{model}.activated-pathways.tsv",
        spia_table_inhibited="results/tables/pathways/{model}.inhibited-pathways.tsv",
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
        "v2.6.0/utils/datavzrd"


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
        "v2.6.0/utils/datavzrd"


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
        "v2.6.0/utils/datavzrd"
