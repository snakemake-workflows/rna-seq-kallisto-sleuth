from pathlib import Path


# topology- and interaction-aware pathway enrichment analysis


# TODO consider cellphonedb for receptor ligand interaction (Sarah Teichmann, Nature Methods?)
rule spia:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        spia_db="resources/spia-db.rds",
    output:
        table="results/tables/pathways/{model}.pathways.tsv",
        table_activated="results/tables/pathways/{model}.activated-pathways.tsv",
        table_inhibited="results/tables/pathways/{model}.inhibited-pathways.tsv",
        plots="results/plots/pathways/{model}.spia-perturbation-plots.pdf",
    params:
        bioc_species_pkg=bioc_species_pkg,
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/tables/pathways/{model}.spia-pathways.log",
    threads: 38
    script:
        "../scripts/spia.R"


## gene set enrichment analysis


rule fgsea:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
    output:
        enrichment=report(
            "results/tables/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        rank_ties=report(
            "results/tables/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        significant=report(
            "results/tables/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        plot=report(
            "results/plots/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        plot_collapsed=report(
            "results/plots/fgsea/{model}.collapsed_pathways.table-plot.pdf",
            caption="../report/fgsea-collapsed-table-plot.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
    params:
        bioc_species_pkg=bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        common_src=str(workflow.source_path("../scripts/common.R")),
    conda:
        enrichment_env
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


rule fgsea_plot_gene_sets:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv",
    output:
        report(
            directory("results/plots/fgsea/{model}"),
            patterns=["{model}.{gene_set}.gene-set-plot.pdf"],
            caption="../report/plot-fgsea-gene-set.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
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
        bioc_species_pkg=bioc_species_pkg,
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
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        enrichment_sig_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
        plot=expand(
            "results/plots/go_terms/{{model}}.go_term_enrichment_{ns}.gene_fdr_{{gene_fdr}}.go_term_fdr_{{go_term_fdr}}.pdf",
            ns=["BP", "CC", "MF"],
        ),
    params:
        species=get_bioc_species_name(),
        model=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        gene_fdr=lambda wc: wc.gene_fdr.replace("-", "."),
        go_term_fdr=lambda wc: wc.go_term_fdr.replace("-", "."),
    conda:
        "../envs/goatools.yaml"
    log:
        "logs/goatools/tables_and_plots.{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
