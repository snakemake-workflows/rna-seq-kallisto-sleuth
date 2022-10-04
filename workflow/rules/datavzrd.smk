rule render_datavzrd_config_spia:
    input:
        template=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        spia_results="results/tables/pathways/{model}.pathways.tsv",
        barplot="results/plots/pathways/{model}.pathways.html",
    output:
        "results/datavzrd/spia/{model}.yaml"
    log:
        "logs/yte/render-datavzrd-config-spia/{model}.log"
    template_engine:
        "yte"

rule datavzrd:
    input:
        config="results/datavzrd/{type}/{model}.yaml",
        # files required for rendering the given configs
        spia_results="results/tables/pathways/{model}.pathways.tsv",
        barplot="results/plots/pathways/{model}.pathways.html",
    output:
        report(
            directory("results/datavzrd-reports/{type}.{model}"),
            htmlindex="index.html",
            # see https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
            # for additional options like caption, categories and labels
        ),
    log:
        "logs/datavzrd-report/{type}.{model}/{type}.{model}.log",
    wrapper:
        "v1.12.2/utils/datavzrd"