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
                "signed_pi_value",
            )
        ]
    )
    return df


prepared_diffexp_x = prepare(diffexp_x)
prepared_diffexp_y = prepare(diffexp_y)
combined = (
    prepared_diffexp_x.join(prepared_diffexp_y, on=["Name"], suffix="_y")
    .with_columns(pl.min_horizontal("Combined FDR", "Combined FDR_y").alias("fdr_min"))
    .with_columns(
        pl.min_horizontal("signed_pi_value", "signed_pi_value_y").alias("min_pi_value")
    )
    .filter(pl.col("fdr_min") <= 0.05)
    .rename(
        {
            "total perturbation accumulation": effect_x,
            "total perturbation accumulation_y": effect_y,
            "fdr_min": "min fdr",
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
print(combined.columns)
combined_sorted = combined.sort(pl.col("min_pi_value").abs(), descending=True)
combined_pd = combined_sorted.select(
    pl.col(
        "Name",
        "min fdr",
        effect_x,
        effect_y,
        "difference",
        "pathway id",
        "min_pi_value",
    )
).to_pandas()
combined_pd.to_csv(snakemake.output[0], sep="\t", index=False)

# we cannot use vegafusion here because it makes the point selection impossible since
# it prunes the required ext_gene column
# alt.data_transformers.enable("vegafusion")
alt.data_transformers.disable_max_rows()


alt.data_transformers.disable_max_rows()
point_selector = alt.selection_point(fields=["Name"], empty=False)

points = (
    alt.Chart(combined_pd)
    .mark_circle(size=15, tooltip={"content": "data"})
    .encode(
        alt.X(
            effect_x,
            scale=alt.Scale(type="symlog", nice=False),
            axis=alt.Axis(grid=False),
        ),
        alt.Y(
            effect_y,
            scale=alt.Scale(type="symlog", nice=False),
            axis=alt.Axis(grid=False),
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

x_axis = (
    alt.Chart(pl.DataFrame({effect_x: [0, 0], effect_y: [min_value, max_value]}))
    .mark_line(color="lightgrey")
    .encode(
        x=effect_x,
        y=effect_y,
        strokeDash=alt.value([5, 5]),
    )
)

y_axis = (
    alt.Chart(pl.DataFrame({effect_x: [min_value, max_value], effect_y: [0, 0]}))
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
        text=alt.condition(point_selector, "Name", alt.value("")),
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
        text=alt.condition(point_selector, "Name", alt.value("")),
    )
)

zero_lines = (
    alt.Chart(pl.DataFrame({"zero": [0]}))
    .mark_rule(color="black")
    .encode(
        x=alt.X("zero", axis=alt.Axis(title="")),
        y=alt.Y("zero", axis=alt.Axis(title="")),
    )
)

chart = (
    alt.layer(x_axis, y_axis, line, points, text_background, text)
    .add_params(point_selector)
    .interactive()
)

chart.save(snakemake.output[1], inline=True)
