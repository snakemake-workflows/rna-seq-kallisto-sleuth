import polars as pl
import altair as alt
import pandas as pd

try:
    input_file = snakemake.input[0]
    df = (
        pl.read_csv(input_file, separator="\t")
        .drop_nulls()
        .filter(~pl.col("pi_value").is_infinite())
    )

    num_go_terms = snakemake.params["num_go_terms"]
    df_top_go_terms = (
        df.sort("pi_value", descending=True)
        .head(num_go_terms)
        .select(["term", "effect", "pi_value"])
    )

    pandas_df = df_top_go_terms.to_pandas()
except Exception as e:
    print(f"File is empty: {e}")
    pandas_df = pd.DataFrame(columns=["term", "effect", "pi_value"])


chart = (
    alt.Chart(pandas_df)
    .mark_bar()
    .encode(
        x=alt.X("effect:Q", title="Effect"),
        y=alt.Y("term:N", sort="-color", title="Term"),
        color=alt.Color(
            "pi_value:Q",
            scale=alt.Scale(scheme="viridis"),
            title="pi_value",
        ),
        tooltip=["term", "effect", "pi_value"],
    )
    .properties(title="Top GO-Terms", width=600, height=400)
)

chart.save(snakemake.output[0])
