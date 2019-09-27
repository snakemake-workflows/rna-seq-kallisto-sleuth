rule spia:
    input:
        diffexp="tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv"
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

## gene set enrichment analysis

checkpoint fgsea:
    input:
        diffexp="tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        gene_sets=config["enrichment"]["gene_sets_file"]
    output:
        enrichment=report(
            "tables/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis"
            ),
        rank_ties=report(
            "tables/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis"
            ),
        significant=report(
            "tables/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis"
            ),
        plot=report(
            "plots/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis"
            )
    params:
        species=config["ref"]["species"],
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        nperm=config["enrichment"]["fgsea"]["nperm"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/fgsea.yaml"
    threads: 8
    script:
        "../scripts/fgsea.R"

rule fgsea_plot_gene_set:
    input:
        diffexp="tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        gene_sets=config["enrichment"]["gene_sets_file"],
        sig_gene_sets="tables/fgsea/{model}.sig-gene-sets.tsv"
    output:
        report(
            "plots/fgsea/{model}.{gene_set}.gene-set-plot.pdf",
            caption="../report/fgsea-gene-set-plot.rst",
            category="Gene set enrichment analysis"
            )
    params:
        model=get_model,
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/fgsea.yaml"
    script:
        "../scripts/fgsea_plot_gene_set.R"

## gene ontology term enrichment analysis

rule biomart_ens_gene_to_go:
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

rule goatools_go_enrichment:
    input:
        obo="data/ref/gene_ontology.obo",
        ens_gene_to_go="data/ref/ens_gene_to_go.tsv",
        diffexp="tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv"
    output:
        enrichment=report(
            "tables/go_terms/{model}.genes-mostsigtrans.diffexp.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
            caption="../report/go-enrichment-mostsigtrans-table.rst",
            category="GO term enrichment analysis"
            ),
        plot=report(
            expand("plots/go_terms/{{model}}.genes-mostsigtrans.diffexp.go_term_enrichment_{ns}.gene_fdr_{{gene_fdr}}.go_term_fdr_{{go_term_fdr}}.pdf",
                    ns = ['BP', 'CC', 'MF']
                    ),
            caption="../report/go-enrichment-mostsigtrans-plot.rst",
            category="GO term enrichment analysis"
            )
    params:
        species=config["ref"]["species"],
        model=get_model,
        gene_fdr=lambda wc: wc.gene_fdr.replace('-','.'),
        go_term_fdr=lambda wc: wc.go_term_fdr.replace('-','.')
    conda:
        "../envs/goatools.yaml"
    script:
        "../scripts/goatools-go-enrichment-analysis.py"



