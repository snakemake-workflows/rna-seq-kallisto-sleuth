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
        #transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        "results/datavzrd/diffexp/{model}.yaml",
    params:
        samples=expand("{unit.sample}-{unit.unit}",unit=units.itertuples()),
    log:
        "logs/yte/render-datavzrd-config-diffexp/{model}.log"
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
        "v1.19.0/utils/datavzrd"

rule diffexp_datavzrd:
    input:
        config="results/datavzrd/diffexp/{model}.yaml",
        #transcripts="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
        genes_aggregated="results/tables/diffexp/{model}.genes-aggregated.diffexp.tsv",
        genes_representative="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        
    output:
        report(
            directory("results/diffexp-reports/diffexp.{model}"),
            htmlindex="index.html",
            # see https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
            # for additional options like caption, categories and labels
        ),
    params:
        samples=expand("{unit.sample}-{unit.unit}",unit=units.itertuples()),
    log:
        "logs/datavzrd-report/diffexp.{model}/diffexp.{model}.log",
    wrapper:
        "v1.19.0/utils/datavzrd"