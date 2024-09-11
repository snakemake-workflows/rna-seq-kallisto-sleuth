from pathlib import Path


# topology- and interaction-aware pathway enrichment analysis


# TODO consider cellphonedb for receptor ligand interaction (Sarah Teichmann, Nature Methods?)
rule spia:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        spia_db="resources/spia-db.rds",
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        table="results/tables/pathways/{model}.pathways.tsv",
        plots="results/plots/pathways/{model}.spia-perturbation-plots.pdf",
    params:
        bioc_species_pkg=bioc_species_pkg,
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    conda:
        enrichment_env
    log:
        "logs/tables/pathways/{model}.spia-pathways.log",
    threads: 38
    script:
        "../scripts/spia.R"
