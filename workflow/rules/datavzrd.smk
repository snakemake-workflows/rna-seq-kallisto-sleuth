rule render_datavzrd_config_spia:
    input:
        template=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        spia_results="results/tables/pathways/{model}.pathways.tsv",
    output:
        "results/datavzrd/spia/{model}.yaml"
    log:
        "logs/yte/render-datavzrd-config-spia/{model}.log"
    template_engine:
        "yte"

rule render_datavzrd_config_diffexp:
    input:
        template=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        "results/datavzrd/diffexp/{model}.yaml",
    params:
        samples=get_model_samples,
    log:
        "logs/yte/render-datavzrd-config-diffexp/{model}.log"
    template_engine:
        "yte"

rule render_datavzrd_config_go_enrichment:
    input:
        template=workflow.source_path("../resources/datavzrd/go-enrichment-template.yaml"),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        sig_go="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        "results/datavzrd/go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.yaml",
    log:
        "logs/yte/render-datavzrd-config-go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log"
    template_engine:
        "yte"

rule spia_datavzrd:
    input:
        config="results/datavzrd/spia/{model}.yaml",
        # files required for rendering the given configs
        spia_results="results/tables/pathways/{model}.pathways.tsv",
        
    output:
        report(
            directory("results/datavzrd-reports/spia-{model}"),
            htmlindex="index.html",
            caption="../report/spia.rst",
            category="Pathway enrichment analysis",
        ),
    log:
        "logs/datavzrd-report/spia-{model}/spia-{model}.log",
    wrapper:
        "v1.19.2/utils/datavzrd"

rule diffexp_datavzrd:
    input:
        config="results/datavzrd/diffexp/{model}.yaml",
        logcount_matrix="results/tables/logcount-matrix/{model}.logcount-matrix.tsv",
        transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        report(
            directory("results/datavzrd-reports/diffexp-{model}"),
            caption="../report/diffexp.rst",
            category="Differential transcript expression",
            htmlindex="index.html",
            # see https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
            # for additional options like caption, categories and labels
        ),
    log:
        "logs/datavzrd-report/diffexp.{model}/diffexp.{model}.log",
    wrapper:
        "v1.19.2/utils/datavzrd"

rule go_enrichment_datavzrd:
    input:
        config="results/datavzrd/go_terms/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.yaml",
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        sig_go="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        report(
            directory("results/datavzrd-reports/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}"),
            htmlindex="index.html",
            caption="../report/go-enrichment-sig_terms.rst",
            category="GO term enrichment analysis",
        ),
    log:
        "logs/datavzrd-report/go_enrichment-{model}/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    wrapper:
        "v1.19.2/utils/datavzrd"