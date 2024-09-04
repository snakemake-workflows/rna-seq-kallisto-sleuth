import pyreadr
import polars as pl
import polars.selectors as cs
import altair as alt

diffexp_x = pl.from_pandas(pyreadr.read_r(snakemake.input[0])[None]).lazy()
diffexp_y = pl.from_pandas(pyreadr.read_r(snakemake.input[1])[None]).lazy()
label_x = list(snakemake.params.labels.keys())[0]
label_y = list(snakemake.params.labels.keys())[1]

effect_x = f"effect {label_x} (beta score)"
effect_y = f"effect {label_y} (beta score)"


def prepare(df):
    return df.select(
        cs.by_name("target_id", "ext_gene", "pval", "qval"),
        cs.starts_with("b_").alias("beta"),
    ).select(
        cs.by_index(range(0, 5))  # only keep first b_ column
    )


combined = (
    prepare(diffexp_x)
    .join(prepare(diffexp_y), on=["target_id", "ext_gene"], suffix="_y")
    .with_columns(
        pl.min_horizontal("qval", "qval_y").alias("qval_min"),
    )
    .filter(pl.col("qval_min") <= 0.05)
    .rename(
        {
            "beta": effect_x,
            "beta_y": effect_y,
            "qval_min": "min q-value",
        }
    )
    .collect()
)

effects = combined.select(pl.col(effect_x, effect_y))
min_value = effects.min().min_horizontal()[0]
max_value = effects.max().max_horizontal()[0]
combined = combined.with_columns(
    abs(pl.col(effect_x) - pl.col(effect_y)).alias("difference")
)
combined_sorted = combined.sort("difference", descending=True)
combined_pd = combined_sorted.select(
    pl.col("ext_gene", "target_id", "min q-value", effect_x, effect_y, "difference")
).to_pandas()
combined_pd.to_csv(snakemake.output[0], sep="\t", index=False)


# we cannot use vegafusion here because it makes the point selection impossible since
# it prunes the required ext_gene column
# alt.data_transformers.enable("vegafusion")
alt.data_transformers.disable_max_rows()

point_selector = alt.selection_point(fields=["ext_gene"], empty=False)

points = (
    alt.Chart(combined_pd)
    .mark_circle(size=15, tooltip={"content": "data"})
    .encode(
        alt.X(effect_x),
        alt.Y(effect_y),
        alt.Color("min q-value", scale=alt.Scale(scheme="viridis")),
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
    alt.Chart(combined_pd)
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
        text=alt.condition(point_selector, "ext_gene", alt.value("")),
    )
)

text = (
    alt.Chart(combined_pd)
    .mark_text(
        align="left",
        baseline="middle",
        dx=5,
        dy=-5,
    )
    .encode(
        x=effect_x,
        y=effect_y,
        text=alt.condition(point_selector, "ext_gene", alt.value("")),
    )
)


chart = (
    alt.layer(line, points, text_background, text)
    .add_params(point_selector)
    .interactive()
)


chart.save(snakemake.output[1], inline=True)
