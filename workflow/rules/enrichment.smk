from pathlib import Path


# topology- and interaction-aware pathway enrichment analysis


# TODO consider cellphonedb for receptor ligand interaction (Sarah Teichmann, Nature Methods?)
rule spia:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        spia_db="resources/spia-db.rds",
    output:
        table=report(
            "results/tables/pathways/{model}.pathways.tsv",
            caption="../report/spia.rst",
            category="Pathway enrichment analysis",
        ),
        plots="results/plots/pathways/{model}.spia-perturbation-plots.pdf",
    params:
        bioc_species_pkg=get_bioc_species_pkg,
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/tables/pathways/{model}.spia-pathways.log",
    threads: 16
    script:
        "../scripts/spia.R"


## gene set enrichment analysis


rule fgsea:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
    output:
        enrichment=report(
            "results/tables/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis",
        ),
        rank_ties=report(
            "results/tables/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis",
        ),
        significant=report(
            "results/tables/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis",
        ),
        plot=report(
            "results/plots/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis",
        ),
        plot_collapsed=report(
            "results/plots/fgsea/{model}.collapsed_pathways.table-plot.pdf",
            caption="../report/fgsea-collapsed-table-plot.rst",
            category="Gene set enrichment analysis",
        ),
    params:
        bioc_species_pkg=get_bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log",
    threads: 8
    script:
        "../scripts/fgsea.R"


rule fgsea_plot_gene_sets:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv",
    output:
        report(
            directory("results/plots/fgsea/{model}"),
            patterns=["{model}.{gene_set}.gene-set-plot.pdf"],
            caption="../report/plot-fgsea-gene-set.rst",
            category="Gene set enrichment analysis",
        ),
    params:
        model=get_model,
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/plots/fgsea/{model}.plot_fgsea_gene_set.log",
    script:
        "../scripts/plot-fgsea-gene-sets.R"


## gene ontology term enrichment analysis


rule ens_gene_to_go:
    output:
        "resources/ontology/ens_gene_to_go.tsv",
    params:
        bioc_species_pkg=get_bioc_species_pkg,
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/resources/ens_gene_to_go.log",
    script:
        "../scripts/ens_gene_to_go.R"


rule download_go_obo:
    output:
        "resources/ontology/gene_ontology.obo",
    params:
        download=config["resources"]["ontology"]["gene_ontology"],
    conda:
        "../envs/curl.yaml"
    log:
        "logs/resources/curl.download_go_obo.log",
    shell:
        "( curl --silent -o {output} {params.download} ) 2> {log}"


rule goatools_go_enrichment:
    input:
        obo="resources/ontology/gene_ontology.obo",
        ens_gene_to_go="resources/ontology/ens_gene_to_go.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
    output:
        enrichment=report(
            "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
            caption="../report/go-enrichment-table.rst",
            category="GO term enrichment analysis",
        ),
        plot=report(
            expand(
                "results/plots/go_terms/{{model}}.go_term_enrichment_{ns}.gene_fdr_{{gene_fdr}}.go_term_fdr_{{go_term_fdr}}.pdf",
                ns=["BP", "CC", "MF"],
            ),
            caption="../report/go-enrichment-plot.rst",
            category="GO term enrichment analysis",
        ),
    params:
        species=get_bioc_species_name(),
        model=get_model,
        gene_fdr=lambda wc: wc.gene_fdr.replace("-", "."),
        go_term_fdr=lambda wc: wc.go_term_fdr.replace("-", "."),
    conda:
        "../envs/goatools.yaml"
    log:
        "logs/goatools/tables_and_plots.{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
