import polars as pl
import polars.selectors as cs
import altair as alt

diffexp_x = pl.read_csv(snakemake.input[0], separator="\t").lazy()
diffexp_y = pl.read_csv(snakemake.input[1], separator="\t").lazy()
label_x = list(snakemake.params.labels[0].keys())[0]
label_y = list(snakemake.params.labels[1].keys())[0]

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
    df = df.select(
        [cs.by_name("GO", "term", "p_uncorrected", "p_fdr_bh", "study_items")]
    )

    # map_elements functions to columns using with_columns
    df = df.with_columns(
        [pl.col("study_items").map_elements(extract_study_items).alias("parsed_terms")]
    )
    df = df.with_columns(
        [
            pl.col("parsed_terms")
            .map_elements(lambda x: calculate_sums(x)[0])
            .alias("cumulative_b_scores_positive"),
            pl.col("parsed_terms")
            .map_elements(lambda x: calculate_sums(x)[1])
            .alias("cumulative_b_scores_negative"),
        ]
    )

    return df


prepared_diffexp_x = prepare(diffexp_x)
prepared_diffexp_y = prepare(diffexp_y)
combined = (
    prepared_diffexp_x.join(prepared_diffexp_y, on=["GO", "term"], suffix="_y")
    .with_columns(pl.min_horizontal("p_fdr_bh", "p_fdr_bh_y").alias("p_fdr_bh_min"))
    .filter(pl.col("p_fdr_bh_min") <= 0.05)
    .rename(
        {
            "cumulative_b_scores_positive": effect_x_pos,
            "cumulative_b_scores_positive_y": effect_y_pos,
            "cumulative_b_scores_negative": effect_x_neg,
            "cumulative_b_scores_negative_y": effect_y_neg,
            "p_fdr_bh_min": "min_p_fdr_bh",
        }
    )
    .collect()
)


# we cannot use vegafusion here because it makes the point selection impossible since
# it prunes the required ext_gene column
# alt.data_transformers.enable("vegafusion")
xlabel = f"effect {label_x}"
ylabel = f"effect {label_y}"
combined = combined.filter(
    (pl.col(effect_x_pos) != 0)
    & (pl.col(effect_y_pos) != 0)
    & (pl.col(effect_x_neg) != 0)
    & (pl.col(effect_y_neg) != 0)
)
combined = combined.with_columns(
    pl.max_horizontal(
        abs(pl.col(effect_x_pos) - pl.col(effect_y_pos)),
        abs(pl.col(effect_x_neg) - pl.col(effect_y_neg)),
    ).alias("difference")
)
combined = combined.with_columns(
    (-pl.col("min_p_fdr_bh").log(base=10) * pl.col("difference")).alias(
        "signed_pi_value"
    )
)
combined_sorted = combined.sort(pl.col("signed_pi_value").abs(), descending=True)
combined_pd = combined_sorted.select(
    pl.col(
        "GO",
        "term",
        "min_p_fdr_bh",
        effect_x_pos,
        effect_y_pos,
        effect_x_neg,
        effect_y_neg,
        "difference",
        "signed_pi_value",
    )
)
combined_pd.to_pandas().to_csv(snakemake.output[0], sep="\t", index=False)

point_selector = alt.selection_point(fields=["term"], empty=False)


def plot(df, effect_x, effect_y, title, xlabel, ylabel):
    # Filter out rows where either effect_x or effect_y is zero because of logarithmic scale
    df = df.select(pl.col("GO", "term", "min_p_fdr_bh", effect_x, effect_y))

    effects = df.select(pl.col(effect_x), pl.col(effect_y))
    min_value = effects.min().min_horizontal()[0]
    max_value = effects.max().max_horizontal()[0]
    df = df.to_pandas()

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
        .properties(title=title)
        .interactive()
    )
    return chart


# Update the plot function calls to include the logarithmic scale and filter out zero values
positive_chart = plot(
    combined_pd, effect_x_pos, effect_y_pos, "Positive effect", xlabel, ylabel
)
negative_chart = plot(
    combined_pd, effect_x_neg, effect_y_neg, "Negative effect", xlabel, ylabel
)

final_chart = alt.hconcat(positive_chart, negative_chart).resolve_scale(color="shared")


final_chart.save(snakemake.output[1], inline=True)
