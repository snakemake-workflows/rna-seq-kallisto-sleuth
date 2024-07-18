rule spia_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
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
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
    wrapper:
        "v3.13.8/utils/datavzrd"


rule diffexp_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/diffexp-template.yaml"),
        # optional files required for rendering the given config
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
    log:
        "logs/datavzrd-report/diffexp.{model}/diffexp.{model}.log",
    params:
        extra="",
        model=get_model,
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        samples=get_model_samples,
    wrapper:
        "v3.13.8/utils/datavzrd"


rule go_enrichment_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/go-enrichment-template.yaml"),
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
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
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        samples=get_model_samples,
    wrapper:
        "v3.13.8/utils/datavzrd"
