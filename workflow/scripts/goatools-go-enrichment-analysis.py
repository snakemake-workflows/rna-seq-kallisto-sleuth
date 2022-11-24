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

#run one time to initialize
go_items = [val for cat in ["BP", "CC", "MF"] for item, val in goeaobj.nsobjgoea[cat].assoc.items()]

#go_file['per'] = go_file.n_genes/go_file.n_go

if goea_results_all:
    goeaobj.wr_tsv(snakemake.output.enrichment, goea_results_all)
else:
    # write empty table to indicate that nothing was found
    with open(snakemake.output.enrichment, "w") as out:
        print(
            "# GO",
            "NS",
            "enrichment",
            "name",
            "ratio_in_study",
            "ratio_in_pop",
            "p_uncorrected",
            "depth",
            "study_count",
            "p_fdr_bh",
            "study_items",
            sep="\t",
            file=out,
        )


# plot results
ensembl_id_to_symbol = dict(zip(all_genes["ens_gene"], all_genes["ext_gene"]))


# from first plot output file name, create generic file name to trigger
# separate plots for each of the gene ontology name spaces
outplot_generic = (
    snakemake.output.plot[0]
    .replace("_BP.", "_{NS}.", 1)
    .replace("_CC.", "_{NS}.", 1)
    .replace("_MF.", "_{NS}.", 1)
)

goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < fdr_level_go_term]

#https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_python.ipynb
if goea_results_sig != "":
    go_sig_terms = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,\
        x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),\
            ], goea_results_sig)), columns = ['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',\
                'n_study', 'n_go'])
    go_sig_terms['gene_ratio'] = go_sig_terms.n_genes/go_sig_terms.n_go
    go_sig_terms.to_csv(snakemake.output.enrichment_sig_terms, sep='\t', index=False)
else:
    no_sig_terms="no significant terms"
    no_sig_terms.to_csv(snakemake.output.enrichment_sig_terms, sep='\t', index=False)
plot_results(
    outplot_generic,
    # use pvals for coloring
    goea_results=goea_results_sig,
    # print general gene symbol instead of Ensembl ID
    id2symbol=ensembl_id_to_symbol,
    # number of genes to print, or True for printing all (including count of genes)
    study_items=True,
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
