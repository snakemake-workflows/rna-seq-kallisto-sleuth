import sys

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")

import pandas as pd
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_results  # , plot_goid2goobj, plot_gos
# read in directed acyclic graph of GO terms / IDs
obodag = GODag(snakemake.input.obo)


# read in mapping gene ids from input to GO terms / IDs
objanno = IdToGosReader(snakemake.input.ens_gene_to_go, godag=obodag)


# extract namespace(?) -> id2gos mapping
ns2assoc = objanno.get_ns2assc()

for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated genes".format(NS=nspc, N=len(id2gos)))

# read gene diffexp table
all_genes = pd.read_table(snakemake.input.diffexp)
# select genes significantly differentially expressed according to BH FDR of sleuth
fdr_level_gene = float(snakemake.params.gene_fdr)
sig_genes = all_genes[all_genes["qval"] < fdr_level_gene]

# initialize GOEA object
fdr_level_go_term = float(snakemake.params.go_term_fdr)
model = snakemake.params.model

goeaobj = GOEnrichmentStudyNS(
    # list of 'population' of genes looked at in total
    pop=all_genes["ens_gene"].tolist(),
    # geneid -> GO ID mapping
    ns2assoc=ns2assoc,
    # ontology DAG
    godag=obodag,
    propagate_counts=False,
    # multiple testing correction method (fdr_bh is false discovery rate control with Benjamini-Hochberg)
    methods=["fdr_bh"],
    # significance cutoff for method named above
    alpha=fdr_level_go_term,
)

goea_results_all = goeaobj.run_study(sig_genes["ens_gene"].tolist())
ensembl_id_to_symbol = dict(zip(all_genes["ens_gene"], all_genes["ext_gene"]))

go_items = [
    val
    for cat in ["BP", "CC", "MF"]
    for item, val in goeaobj.ns2objgoea[cat].assoc.items()
]

# goea_results_all.rename(columns={'# GO': 'GO'})


if goea_results_all:
    go_terms = pd.DataFrame(
        list(
            map(
                lambda x: [
                    x.GO,
                    x.goterm.name,
                    x.goterm.namespace,
                    x.p_uncorrected,
                    x.p_fdr_bh,
                    x.ratio_in_study,
                    x.ratio_in_pop,
                    x.depth,
                    x.study_count,
                    x.study_items,
                ],
                goea_results_all,
            )
        ),
        columns=[
            "GO",
            "term",
            "class",
            "p_uncorrected",
            "p_fdr_bh",
            "ratio_in_study",
            "ratio_in_pop",
            "depth",
            "study_count",
            "study_items",
        ],
    )
    go_terms["study_items"] = go_terms["study_items"].str.join(",")
    go_terms["study_items"] = go_terms.study_items.str.replace(
        "\w+(?=,|$)", lambda m: ensembl_id_to_symbol.get(m.group(0))
    )
    go_terms.to_csv(snakemake.output.enrichment, sep="\t", index=False)
else:
    # write empty table to indicate that nothing was found
    with open(snakemake.output.enrichment, "w") as out:
        print(
            "GO",
            "term",
            "class",
            "name",
            "p_uncorrected",
            "p_fdr_bh",
            "ratio_in_study",
            "ratio_in_pop",
            "depth",
            "study_count",
            "study_items",
            sep="\t",
            file=out,
        )


# plot results

# from first plot output file name, create generic file name to trigger
# separate plots for each of the gene ontology name spaces
outplot_generic = (
    snakemake.output.plot[0]
    .replace("_BP.", "_{NS}.", 1)
    .replace("_CC.", "_{NS}.", 1)
    .replace("_MF.", "_{NS}.", 1)
)

goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < fdr_level_go_term]

# https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_python.ipynb
if goea_results_sig:
    go_sig_terms = pd.DataFrame(
        list(
            map(
                lambda x: [
                    x.GO,
                    x.goterm.name,
                    x.goterm.namespace,
                    x.p_uncorrected,
                    x.p_fdr_bh,
                    x.ratio_in_study[0],
                    x.ratio_in_study[1],
                    x.ratio_in_pop[0],
                    x.study_items,
                ],
                goea_results_sig,
            )
        ),
        columns=[
            "GO",
            "term",
            "class",
            "p-value",
            "p_corrected",
            "n_genes_diff_exp",
            "n_genes_diff_exp_study",
            "n_go_genes",
            "study_items",
        ],
    )
    go_sig_terms["gene_ratio"] = go_sig_terms.n_genes_diff_exp / go_sig_terms.n_go_genes
    go_sig_terms_sorted = go_sig_terms.sort_values(by=["class", "p_corrected"])
    go_sig_terms_sorted["study_items"] = go_sig_terms_sorted["study_items"].str.join(
        ","
    )
    go_sig_terms_sorted["study_items"] = go_sig_terms_sorted.study_items.str.replace(
        "\w+(?=,|$)", lambda m: ensembl_id_to_symbol.get(m.group(0))
    )
    # Append fold change values to gene names
    gene_to_fold_change = dict(
        zip(sig_genes["ext_gene"], sig_genes.filter(regex=("b_" + model)).iloc[:, 0])
    )
    go_sig_terms_sorted["study_items"] = go_sig_terms_sorted.study_items.astype("str")
    go_sig_terms_sorted["study_items"] = go_sig_terms_sorted.study_items.str.split(",")
    # As go terms contains 0 genes associated with terms, assign empty key to dict fold change as 0
    gene_to_fold_change[""] = 0
    gene_to_fold_change["nan"] = 0
    go_sig_terms_sorted["study_items"] = go_sig_terms_sorted.study_items.map(
        lambda x: ", ".join(
            [f"{gene_symbol}:{gene_to_fold_change[gene_symbol]}" for gene_symbol in x]
        )
    )
    go_sig_terms_sorted.to_csv(
        snakemake.output.enrichment_sig_terms, sep="\t", index=False
    )
else:
    # write empty table to indicate that nothing was found
    with open(snakemake.output.enrichment_sig_terms, "w") as out:
        print(
            "GO",
            "term",
            "class",
            "p-value",
            "p_corrected",
            "n_genes_diff_exp",
            "n_genes_diff_exp_study",
            "n_go_genes",
            "gene_ratio",
            "study_items",
            sep="\t",
            file=out,
        )
plot_results(
    outplot_generic,
    # use pvals for coloring
    goea_results=goea_results_sig,
    # print general gene symbol instead of Ensembl ID
    id2symbol=ensembl_id_to_symbol,
    # number of genes to print, or True for printing all (including count of genes)
    study_items=None,
    # number of genes to print per line
    items_p_line=6,
    # p-values to use for coloring of GO term nodes (argument name determined from code and testing against value "p_uncorrected")
    pval_name="p_fdr_bh",
)

# for all name spaces
for ns in ns2assoc.keys():
    # check if no GO terms were found to be significant
    if len([r for r in goea_results_sig if r.NS == ns]) == 0:
        fig = plt.figure(figsize=(12, 2))
        text = fig.text(
            0.5,
            0.5,
            "No plot generated, because no GO terms were found significant\n"
            "for name space {} and significance levels: genes ({}), GO terms ({}).\n"
            "You might want to check those levels and/or your intermediate data.".format(
                ns, fdr_level_gene, fdr_level_go_term
            ),
            ha="center",
            va="center",
            size=20,
        )
        fig.savefig(outplot_generic.replace("_{NS}.", "_{}.".format(ns)))
