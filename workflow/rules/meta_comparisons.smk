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
<<<<<<< Updated upstream
    notebook:
        "../scripts/compare_diffexp.py.ipynb"
=======
    script:
        "../scripts/compare_diffexp.py"


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
>>>>>>> Stashed changes
