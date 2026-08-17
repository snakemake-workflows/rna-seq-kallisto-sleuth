from pathlib import Path

# topology- and interaction-aware pathway enrichment analysis


# TODO consider cellphonedb for receptor ligand interaction (Sarah Teichmann, Nature Methods?)
rule spia:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        spia_db="resources/spia-db.{database}.rds",
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        table="results/tables/pathways/{model}.{database}.pathways.tsv",
        plots="results/plots/pathways/{model}.{database}.spia-perturbation-plots.pdf",
    log:
        "logs/tables/pathways/{model}.{database}.spia-pathways.log",
    conda:
        enrichment_env
    threads: 32
    resources:
        mem_mb=32000,
    params:
        bioc_species_pkg=bioc_species_pkg,
        covariate=lookup(within=config, dpath="diffexp/models/{model}/primary_variable"),
    script:
        "../scripts/spia.R"


## gene set enrichment analysis


rule fgsea:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        common_src=workflow.source_path("../scripts/common.R"),
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
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log",
    conda:
        enrichment_env
    threads: 25
    params:
        bioc_species_pkg=bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    script:
        "../scripts/fgsea.R"


rule fgsea_plot_gene_sets:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv",
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        report(
            directory("results/plots/fgsea/{model}"),
            patterns=["{model}.{gene_set}.gene-set-plot.pdf"],
            caption="../report/plot-fgsea-gene-set.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
    log:
        "logs/plots/fgsea/{model}.plot_fgsea_gene_set.log",
    conda:
        enrichment_env
    params:
        model=get_model,
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    script:
        "../scripts/plot-fgsea-gene-sets.R"


## gene ontology term enrichment analysis


rule ens_gene_to_go:
    input:
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        "resources/ontology/ens_gene_to_go.tsv",
    log:
        "logs/resources/ens_gene_to_go.log",
    conda:
        enrichment_env
    params:
        bioc_species_pkg=bioc_species_pkg,
    script:
        "../scripts/ens_gene_to_go.R"


rule download_go_obo:
    output:
        "resources/ontology/gene_ontology.obo",
    log:
        "logs/resources/curl.download_go_obo.log",
    conda:
        "../envs/curl.yaml"
    params:
        download=config["resources"]["ontology"]["gene_ontology"],
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
    log:
        "logs/goatools/tables_and_plots.{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    conda:
        "../envs/goatools.yaml"
    params:
        species=get_bioc_species_name(),
        model=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        gene_fdr=lambda wc: wc.gene_fdr.replace("-", "."),
        go_term_fdr=lambda wc: wc.go_term_fdr.replace("-", "."),
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
