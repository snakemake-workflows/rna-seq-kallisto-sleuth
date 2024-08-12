rule meta_compare_diffexp:
    input:
        expand(
            "results/sleuth/diffexp/{model}.genes-representative.diffexp.rds",
            model=lookup(
                dpath="meta_comparisons/comparisons/{meta_comp}/items/*/*",
                within=config,
            ),
        ),
    output:
        "results/tables/diffexp/meta_compare_{meta_comp}.tsv",
        "results/meta_comparison/diffexp/{meta_comp}.json",
    log:
        notebook="logs/meta_compare_diffexp/{meta_comp}.ipynb",
    params:
        labels=lookup(
            dpath="meta_comparisons/comparisons/{meta_comp}/items/*",
            within=config,
        ),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/compare_diffexp.py"


rule meta_compare_enrichment:
    input:
        expand(
            "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
            model=lookup(
                dpath="meta_comparisons/comparisons/{meta_comp}/items/*/*",
                within=config,
            ),
            gene_fdr=str(config["enrichment"]["goatools"]["fdr_genes"]).replace(
            ".", "-"
            ),
            go_term_fdr=str(config["enrichment"]["goatools"]["fdr_go_terms"]).replace(
                ".", "-"
            ),
        ),
    output:
        "results/tables/go_terms/meta_compare_{meta_comp}.tsv",
        "results/meta_comparison/go_terms/{meta_comp}.json",
    log:
        notebook="logs/meta_compare_enrichment/{meta_comp}.ipynb",
    params:
        labels=lookup(
            dpath="meta_comparisons/comparisons/{meta_comp}/items/*",
            within=config,
        ),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/compare_enrichment.py"


rule meta_compare_pathways:
    input:
        expand(
            "results/tables/pathways/{model}.pathways.tsv",
            model=lookup(
                dpath="meta_comparisons/comparisons/{meta_comp}/items/*/*",
                within=config,
            ),
        ),
    output:
        "results/tables/pathways/meta_compare_{meta_comp}.tsv",
        "results/meta_comparison/pathways/{meta_comp}.json",
    log:
        notebook="logs/meta_compare_pathways/{meta_comp}.ipynb",
    params:
        labels=lookup(
            dpath="meta_comparisons/comparisons/{meta_comp}/items/*",
            within=config,
        ),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/compare_pathways.py"
