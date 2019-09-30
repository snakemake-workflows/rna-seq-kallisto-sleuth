include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/quant.smk"
include: "rules/diffexp.smk"
include: "rules/enrichment.smk"

rule all:
    input:
        expand(
            [
                "tables/diffexp/{model}.transcripts.diffexp.tsv",
                "plots/diffexp-heatmap/{model}.diffexp-heatmap.pdf",
                "tables/tpm-matrix/{model}.tpm-matrix.tsv",
                "tables/pathways/{model}.pathways.tsv"
            ],
            model=config["diffexp"]["models"]
        ),
        # results goatools
        expand(
            [
                "tables/go_terms/{model}.genes-mostsigtrans.diffexp.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
                "plots/go_terms/{model}.genes-mostsigtrans.diffexp.go_term_enrichment_{go_ns}.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.pdf"
            ],
            model=config["diffexp"]["models"],
            go_ns=["BP", "CC", "MF"],
            gene_fdr=str(config["enrichment"]["goatools"]["fdr_genes"]).replace('.','-'),
            go_term_fdr=str(config["enrichment"]["goatools"]["fdr_go_terms"]).replace('.','-')
        ),
        # results fgsea
        expand(
            [
                "tables/fgsea/{model}.all-gene-sets.tsv",
                "tables/fgsea/{model}.sig-gene-sets.tsv",
                "plots/fgsea/{model}.table-plot.pdf"
            ],
            model=config["diffexp"]["models"],
            gene_set_fdr=str(config["enrichment"]["fgsea"]["fdr_gene_set"]).replace('.','-'),
            nperm=str(config["enrichment"]["fgsea"]["nperm"])
        ),
        [get_fgsea_plots(model) for model in config["diffexp"]["models"]],
        # sleuth p-value histogram plots
        expand("plots/diffexp/{model}.{level}.diffexp-pval-hist.pdf",
                model=config["diffexp"]["models"],
                level=["transcripts", "genes-aggregated", "genes-mostsigtrans" ]
        ),
        # PCA plots of kallisto results, each coloured for a different covariate
        expand("plots/pca/{covariate}.pca.pdf", covariate=samples.columns[samples.columns != "sample"]),
        [get_bootstrap_plots(model) for model in config["diffexp"]["models"]],
        [get_bootstrap_plots(model, config["bootstrap_plots"]["genes_of_interest"])
            for model in config["diffexp"]["models"] ],
        expand("plots/fld/{unit.sample}-{unit.unit}.fragment-length-dist.pdf", unit=units[["sample", "unit"]].itertuples())
