import polars as pl
import polars.selectors as cs
import altair as alt


diffexp_x = pl.read_csv(snakemake.input[0], separator="\t").lazy()
diffexp_y = pl.read_csv(snakemake.input[1], separator="\t").lazy()
label_x = list(snakemake.params.labels[0].keys())[0]
label_y = list(snakemake.params.labels[1].keys())[0]

effect_x = f"effect {label_x}"
effect_y = f"effect {label_y}"


def prepare(df):
    # Select necessary columns and filter
    df = df.select(
        [
            cs.by_name(
                "Name",
                "Combined p-value",
                "Combined FDR",
                "total perturbation accumulation",
                "pathway id",
            )
        ]
    )
    return df


prepared_diffexp_x = prepare(diffexp_x)
prepared_diffexp_y = prepare(diffexp_y)
combined = (
    prepared_diffexp_x.join(prepared_diffexp_y, on=["Name"], suffix="_y")
    .with_columns(
        pl.col("Combined FDR").cast(pl.Float64),
        pl.col("Combined FDR_y").cast(pl.Float64),
        pl.col("total perturbation accumulation").cast(pl.Float64),
        pl.col("total perturbation accumulation_y").cast(pl.Float64),
    )
    .with_columns(pl.min_horizontal("Combined FDR", "Combined FDR_y").alias("min fdr"))
    .filter(pl.col("min fdr") <= 0.05)
    .rename(
        {
            "total perturbation accumulation": effect_x,
            "total perturbation accumulation_y": effect_y,
        }
    )
    .collect()
)

if not combined.is_empty():
    combined = (
        combined.with_columns(
            abs(pl.col(effect_x) - pl.col(effect_y)).alias("difference")
        )
        .with_columns(
            (-pl.col("min fdr").log(base=10) * pl.col("difference")).alias("pi_value")
        )
        .sort(pl.col("pi_value").abs(), descending=True)
    )
else:
    combined = combined.with_columns(
        [pl.lit(None).alias("difference"), pl.lit(None).alias("pi_value")]
    )
combined_pd = combined.select(
    pl.col(
        "Name",
        "min fdr",
        effect_x,
        effect_y,
        "difference",
        "pathway id",
        "pi_value",
    )
)
combined_pd.to_pandas().to_csv(snakemake.output[0], sep="\t", index=False)

df = combined_pd.select(pl.col("min fdr", effect_x, effect_y)).to_pandas()
min_value = min(df[effect_x].min(), df[effect_y].min())
max_value = max(df[effect_x].max(), df[effect_y].max())
point_selector = alt.selection_single(fields=["term"], empty="all")

alt.data_transformers.disable_max_rows()
points = (
    alt.Chart(df)
    .mark_circle(size=15, tooltip={"content": "data"})
    .encode(
        alt.X(
            effect_x,
            title=label_x,
            scale=alt.Scale(type="symlog", nice=False),
            axis=alt.Axis(grid=True),
        ),
        alt.Y(
            effect_y,
            title=label_y,
            scale=alt.Scale(type="symlog", nice=False),
            axis=alt.Axis(grid=True),
        ),
        alt.Color("min fdr", scale=alt.Scale(scheme="viridis")),
        opacity=alt.value(0.5),
    )
)

line = (
    alt.Chart(
        pl.DataFrame(
            {effect_x: [min_value, max_value], effect_y: [min_value, max_value]}
        )
    )
    .mark_line(color="lightgrey")
    .encode(
        x=effect_x,
        y=effect_y,
        strokeDash=alt.value([5, 5]),
    )
)

text_background = (
    alt.Chart(df)
    .mark_text(
        align="left",
        baseline="middle",
        dx=5,
        dy=-5,
        fill="white",
        stroke="white",
        strokeWidth=5,
    )
    .encode(
        x=effect_x,
        y=effect_y,
        text=alt.condition(point_selector, "term", alt.value("")),
    )
)

text = (
    alt.Chart(df)
    .mark_text(
        align="left",
        baseline="middle",
        dx=5,
        dy=-5,
    )
    .encode(
        x=effect_x,
        y=effect_y,
        text=alt.condition(point_selector, "term", alt.value("")),
    )
)

chart = (
    alt.layer(line, points, text_background, text)
    .add_params(point_selector)
    .interactive()
)

chart.save(snakemake.output[1], inline=True)
