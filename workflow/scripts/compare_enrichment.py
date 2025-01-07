import polars as pl
import polars.selectors as cs
import altair as alt
import sys

sys.stderr = open(snakemake.log[0], "w")

diffexp_x = pl.read_csv(snakemake.input[0], separator="\t").lazy()
diffexp_y = pl.read_csv(snakemake.input[1], separator="\t").lazy()
label_x = list(snakemake.params.labels)[0]
label_y = list(snakemake.params.labels)[1]

effect_x_pos = f"positive effect {label_x}"
effect_y_pos = f"positive effect {label_y}"
effect_x_neg = f"negative effect {label_x}"
effect_y_neg = f"negative effect {label_y}"


def extract_study_items(value):
    if value and value.strip() != "":
        gene_values = value.split(", ")
        data = []
        for gene_value in gene_values:
            parts = gene_value.split(":")
            gene = parts[0]
            val = float(parts[1])
            data.append({"gene": gene, "value": val})
        return data
    else:
        return []


def calculate_sums(parsed_terms):
    positive_sum = sum(item["value"] for item in parsed_terms if item["value"] > 0)
    negative_sum = sum(item["value"] for item in parsed_terms if item["value"] < 0)
    return positive_sum, negative_sum


def prepare(df):
    # Select necessary columns
    df = (
        df.select(
            [cs.by_name("GO", "term", "p_uncorrected", "p_fdr_bh", "study_items")]
        )
        .with_columns(
            [
                pl.col("study_items")
                .map_elements(
                    extract_study_items,
                    return_dtype=pl.List(
                        pl.Struct(
                            [pl.Field("gene", pl.Utf8), pl.Field("value", pl.Float64)]
                        )
                    ),
                )
                .alias("parsed_terms")
            ]
        )
        .with_columns(
            [
                pl.col("parsed_terms")
                .map_elements(lambda x: calculate_sums(x)[0], return_dtype=pl.Float64)
                .alias("cumulative_b_scores_positive"),
                pl.col("parsed_terms")
                .map_elements(lambda x: calculate_sums(x)[1], return_dtype=pl.Float64)
                .alias("cumulative_b_scores_negative"),
            ]
        )
    )

    return df


prepared_diffexp_x = prepare(diffexp_x)
prepared_diffexp_y = prepare(diffexp_y)
combined = (
    prepared_diffexp_x.join(
        prepared_diffexp_y, on=["GO", "term"], how="outer", suffix="_y"
    )
    .rename(
        {
            "cumulative_b_scores_positive": effect_x_pos,
            "cumulative_b_scores_positive_y": effect_y_pos,
            "cumulative_b_scores_negative": effect_x_neg,
            "cumulative_b_scores_negative_y": effect_y_neg,
        }
    )
    .with_columns(
        pl.col(effect_x_pos).cast(pl.Float64),
        pl.col(effect_y_pos).cast(pl.Float64),
        pl.col(effect_x_neg).cast(pl.Float64),
        pl.col(effect_y_neg).cast(pl.Float64),
        pl.col("p_fdr_bh").cast(pl.Float64),
        pl.col("p_fdr_bh_y").cast(pl.Float64),
    )
    .with_columns(pl.min_horizontal("p_fdr_bh", "p_fdr_bh_y").alias("min_p_fdr_bh"))
    .fill_null(0)  # Set missing values to 0
    .filter(pl.col("min_p_fdr_bh") <= 0.05)
    .collect()
)


# we cannot use vegafusion here because it makes the point selection impossible since
# it prunes the required ext_gene column
# alt.data_transformers.enable("vegafusion")


if not combined.is_empty():
    combined = (
        combined.with_columns(
            pl.max_horizontal(
                abs(pl.col(effect_x_pos) - pl.col(effect_y_pos)),
                abs(pl.col(effect_x_neg) - pl.col(effect_y_neg)),
            ).alias("difference")
        )
        .with_columns(
            (-pl.col("min_p_fdr_bh").log(base=10) * pl.col("difference")).alias(
                "pi_value"
            )
        )
        .sort(pl.col("pi_value").abs(), descending=True)
    )
else:
    combined = combined.with_columns(
        [pl.lit(None).alias("difference"), pl.lit(None).alias("pi_value")]
    )


combined_pd = combined.select(
    pl.col(
        "GO",
        "term",
        "min_p_fdr_bh",
        effect_x_pos,
        effect_y_pos,
        effect_x_neg,
        effect_y_neg,
        "difference",
        "pi_value",
    )
).to_pandas()

combined_pd.to_csv(snakemake.output[0], sep="\t", index=False)


def plot(df, effect_x, effect_y, title, xlabel, ylabel):
    # Filter out rows where either effect_x or effect_y is zero because of logarithmic scale
    min_value = min(df[effect_x].min(), df[effect_y].min())
    max_value = max(df[effect_x].max(), df[effect_y].max())
    point_selector = alt.selection_single(fields=["term"], empty=False)

    alt.data_transformers.disable_max_rows()
    points = (
        alt.Chart(df)
        .mark_circle(size=15, tooltip={"content": "data"})
        .encode(
            alt.X(
                effect_x,
                title=xlabel,
                scale=alt.Scale(type="symlog", nice=False),
                axis=alt.Axis(grid=True),
            ),
            alt.Y(
                effect_y,
                title=ylabel,
                scale=alt.Scale(type="symlog", nice=False),
                axis=alt.Axis(grid=True),
            ),
            alt.Color("min_p_fdr_bh", scale=alt.Scale(scheme="viridis")),
            opacity=alt.value(0.5),
        )
    )

    min_value = min(df[effect_x].min(), df[effect_y].min())
    max_value = max(df[effect_x].max(), df[effect_y].max())
    line = (
        alt.Chart(
            pl.DataFrame(
                {effect_x: [min_value, max_value], effect_y: [min_value, max_value]},
                schema={effect_x: pl.Float64, effect_y: pl.Float64},
            )
        )
        .mark_line(color="lightgrey")
        .encode(
            x=effect_x,
            y=effect_y,
            strokeDash=alt.value([5, 5]),
        )
    )

    point_selector = alt.selection_single(fields=["term"], empty=False)
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
        # alt.layer(line, points, text_background, text)
        alt.layer(line, points, text_background, text)
        .add_params(point_selector)
        .properties(title=title)
        .interactive()
    )
    return chart


prepared_diffexp_x = prepare(diffexp_x)
prepared_diffexp_y = prepare(diffexp_y)
combined = (
    prepared_diffexp_x.join(prepared_diffexp_y, on=["GO", "term"], suffix="_y")
    .rename(
        {
            "cumulative_b_scores_positive": effect_x_pos,
            "cumulative_b_scores_positive_y": effect_y_pos,
            "cumulative_b_scores_negative": effect_x_neg,
            "cumulative_b_scores_negative_y": effect_y_neg,
        }
    )
    .cast(
        {
            cs.by_name(
                effect_x_pos,
                effect_y_pos,
                effect_x_neg,
                effect_y_neg,
                "p_fdr_bh",
                "p_fdr_bh_y",
            ): pl.Float64
        }
    )
    .with_columns(pl.min_horizontal("p_fdr_bh", "p_fdr_bh_y").alias("min_p_fdr_bh"))
    .filter(pl.col("min_p_fdr_bh") <= 0.05)
    .filter(
        (pl.col(effect_x_pos) != 0)
        & (pl.col(effect_y_pos) != 0)
        & (pl.col(effect_x_neg) != 0)
        & (pl.col(effect_y_neg) != 0)
    )
    .collect()
)


# we cannot use vegafusion here because it makes the point selection impossible since
# it prunes the required ext_gene column
# alt.data_transformers.enable("vegafusion")


if not combined.is_empty():
    combined = (
        combined.with_columns(
            pl.max_horizontal(
                abs(pl.col(effect_x_pos) - pl.col(effect_y_pos)),
                abs(pl.col(effect_x_neg) - pl.col(effect_y_neg)),
            ).alias("difference")
        )
        .with_columns(
            (-pl.col("min_p_fdr_bh").log(base=10) * pl.col("difference")).alias(
                "pi_value"
            )
        )
        .sort(pl.col("pi_value").abs(), descending=True)
    )
else:
    combined = combined.with_columns(
        [pl.lit(None).alias("difference"), pl.lit(None).alias("pi_value")]
    )


combined = combined.select(
    [
        "GO",
        "term",
        "min_p_fdr_bh",
        effect_x_pos,
        effect_y_pos,
        effect_x_neg,
        effect_y_neg,
        "difference",
        "pi_value",
    ]
)

combined.write_csv(snakemake.output[0], separator="\t")


# Update the plot function calls to include the logarithmic scale and filter out zero values
positive_chart = plot(
    combined_pd,
    effect_x_pos,
    effect_y_pos,
    "Positive effect",
    f"effect {label_x}",
    f"effect {label_y}",
)
negative_chart = plot(
    combined_pd,
    effect_x_neg,
    effect_y_neg,
    "Negative effect",
    f"effect {label_x}",
    f"effect {label_y}",
)

final_chart = alt.hconcat(positive_chart, negative_chart).resolve_scale(color="shared")


final_chart.save(snakemake.output[1], inline=True)
