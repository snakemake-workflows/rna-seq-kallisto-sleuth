{{ snakemake.params.desc }}
Each point represents a gene, x axis shows the {{ snakemake.params.labels[0] }} effect, y-axis the {{ snakemake.params.labels[1] }} effect (both as log2 fold change).
The color encodes the corresponding q-value.
By clicking on points, their label can be displayed.
Holding the Shift key allows to select or deselect labels for multiple genes.