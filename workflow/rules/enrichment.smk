
rule download_bioconductor_species_database:
    output:
        "resources/bioconductor/{package}/DESCRIPTION"
    params:
        lib = lambda wc, output: output[0].split(sep=wc.package+"/DESCRIPTION")[0]
    conda:
        "../envs/biocmanager.yaml"
    log:
        "logs/resources/bioconductor/{package}.log"
    script:
        "../scripts/download_bioconductor_package.R"


# topology- and interaction-aware pathway enrichment analysis

rule spia:
    input:
        samples="results/sleuth/samples.tsv",
        species_anno=get_bioc_pkg_path,
        diffexp="results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv"
    output:
        "results/tables/pathways/{model}.pathways.tsv"
    params:
        bioc_pkg=get_bioc_species_pkg,
        species=config["resources"]["ref"]["species"],
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/spia.yaml"
    log:
        "logs/tables/pathways/{model}.spia-pathways.log"
    threads: 16
    script:
        "../scripts/spia.R"

## gene set enrichment analysis

checkpoint fgsea:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        species_anno=get_bioc_pkg_path,
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"]
    output:
        enrichment=report(
            "results/tables/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis"
        ),
        rank_ties=report(
            "results/tables/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis"
        ),
        significant=report(
            "results/tables/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis"
        ),
        plot=report(
            "results/plots/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis"
        )
    params:
        bioc_pkg=get_bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        nperm=config["enrichment"]["fgsea"]["nperm"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/fgsea.yaml"
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log"
    threads: 8
    script:
        "../scripts/fgsea.R"

rule fgsea_plot_gene_set:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv"
    output:
        report(
            "results/plots/fgsea/{model}.{gene_set}.gene-set-plot.pdf",
            caption="../report/fgsea-gene-set-plot.rst",
            category="Gene set enrichment analysis"
        )
    params:
        model=get_model,
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/fgsea.yaml"
    log:
        "logs/plots/fgsea/{model}.{gene_set}.plot_gene_set.log"
    script:
        "../scripts/fgsea_plot_gene_set.R"

## gene ontology term enrichment analysis

rule ens_gene_to_go:
    input:
        species_anno=get_bioc_pkg_path
    output:
        "resources/ontology/ens_gene_to_go.tsv"
    params:
        bioc_pkg=get_bioc_species_pkg
    conda:
        "../envs/ens_gene_to_go.yaml"
    log:
        "logs/resources/ens_gene_to_go.log"
    script:
        "../scripts/ens_gene_to_go.R"


rule download_go_obo:
    output:
        "resources/ontology/gene_ontology.obo"
    params:
        download=config["resources"]["ontology"]["gene_ontology"]
    conda:
        "../envs/curl.yaml"
    log:
        "logs/resources/curl.download_go_obo.log"
    shell:
        "( curl --silent -o {output} {params.download} ) 2> {log}"

rule goatools_go_enrichment:
    input:
        obo="resources/ontology/gene_ontology.obo",
        ens_gene_to_go="resources/ontology/ens_gene_to_go.tsv",
        diffexp="results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv"
    output:
        enrichment=report(
            "results/tables/go_terms/{model}.genes-mostsigtrans.diffexp.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
            caption="../report/go-enrichment-mostsigtrans-table.rst",
            category="GO term enrichment analysis"
        ),
        plot=report(
            expand("results/plots/go_terms/{{model}}.genes-mostsigtrans.diffexp.go_term_enrichment_{ns}.gene_fdr_{{gene_fdr}}.go_term_fdr_{{go_term_fdr}}.pdf",
                    ns = ['BP', 'CC', 'MF']
            ),
            caption="../report/go-enrichment-mostsigtrans-plot.rst",
            category="GO term enrichment analysis"
        )
    params:
        species=config["resources"]["ref"]["species"],
        model=get_model,
        gene_fdr=lambda wc: wc.gene_fdr.replace('-','.'),
        go_term_fdr=lambda wc: wc.go_term_fdr.replace('-','.')
    conda:
        "../envs/goatools.yaml"
    log:
        "logs/goatools/tables_and_plots.{model}.genes-mostsigtrans.diffexp.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.log"
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
