rule spia:
    input:
        diffexp="tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        sleuth="sleuth/{model}.rds"
    output:
        "tables/pathways/{model}.pathways.tsv"
    params:
        species=config["ref"]["species"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/spia.yaml"
    threads: 16
    script:
        "../scripts/spia.R"


rule biomart-ens_gene_to_go:
    output:
        "data/ref/ens_gene_to_go.tsv"
    params:
        species=config["ref"]["species"]
    conda:
        "../envs/biomart-download.yaml"
    script:
        "../scripts/biomart-ens_gene_to_go.R"


rule download_go_obo:
    output:
        "data/ref/gene_ontology.obo"
    params:
        download=config["ref"]["gene_ontology"]
    conda:
        "../envs/curl.yaml"
    log:
        "logs/curl/download_go_obo.log"
    shell:
        "( curl --silent -o {output} {params.download} ) 2> {log}"

