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
