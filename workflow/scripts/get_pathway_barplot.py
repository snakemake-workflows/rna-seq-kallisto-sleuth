import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

pathway_file = pd.read_csv(snakemake.input["pathway_file"], sep="\t")

barplot = alt.Chart(pathway_file).mark_bar().encode(
    x=alt.X("observed total perturbation accumulation:Q"),
    y=alt.Y("Name", sort='-x'), color='Status:N').properties(
        width=600, title = "pathways observed total perturbation accumulation from SPIA of False Discovery Rate global p-values <=0.05" ).transform_filter(
    alt.FieldLTEPredicate(field='False Discovery Rate global p-values', lte = 0.05)
)

barplot.save(snakemake.output[0])
