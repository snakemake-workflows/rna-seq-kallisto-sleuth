import pandas as pd
import altair as alt


def plot(df_plot, effect_x, effect_y, title, selector, color_scheme):
    alt.data_transformers.disable_max_rows()
    points = (
        alt.Chart(df_plot, title=title)
        .mark_circle(size=60)
        .encode(
            alt.Y(
                effect_y,
                title=effect_y,
                scale=alt.Scale(type="log"),
                axis=alt.Axis(grid=False),
            ),
            alt.X(
                effect_x,
                title=effect_x,
                scale=alt.Scale(type="log"),
                axis=alt.Axis(grid=False),
            ),
            color=alt.condition(
                selector,
                alt.Color(
                    f"{name}:N", scale=alt.Scale(scheme=color_scheme), legend=None
                ),
                alt.value("lightgray"),
            ),
            tooltip=[f"{name}:N", f"{effect_x}:Q", f"{effect_y}:Q"],
        )
        .add_params(selector)
    )
    return points


def plot_legend(df_plot, name, selector):
    legend = (
        alt.Chart(df_plot)
        .mark_circle(size=60)
        .encode(
            alt.Y(
                f"{name}",
                axis=alt.Axis(
                    title="", ticks=False, domain=False, labelLimit=700, orient="right"
                ),
            ),
            color=f"{name}",
        )
        .transform_filter(selector)
    )
    return legend


color_scheme = [
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
    "#b15928",
    "#a6cee3",
    "#b2df8a",
    "#fb9a99",
    "#fdbf6f",
    "#cab2d6",
    "#ffff33",
    "#8dd3c7",
    "#ffffb3",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd",
    "#ccebc5",
    "#ffed6f",
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494",
    "#b3b3b3",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ffcc00",
    "#d73027",
    "#4575b4",
    "#91cf60",
    "#e08214",
    "#7fc97f",
    "#beaed4",
    "#fdc086",
    "#ffff99",
    "#386cb0",
    "#f0027f",
    "#bf5b17",
    "#666666",
    "#1b9e77",
    "#d95f02",
]

df = pd.read_csv(snakemake.input[0], sep="\t")
df = df.head(50)
name = snakemake.params["name"]
effect_x = snakemake.params["effect_x"]
effect_y = snakemake.params["effect_y"]

df_positive = df[df[f"{effect_x}"] >= 0]
df_negative = df[df[f"{effect_x}"] < 0]
df_negative[f"{effect_x}"] = df_negative[f"{effect_x}"].abs()

point_selector = alt.selection_point(empty=False)

if df_negative.empty:
    scatter = plot(df_positive, effect_x, effect_y, "", point_selector, color_scheme)
    # Important: You need to copy the df in order to have different datasets, else vega does not bin the plot to a dataset
    legend = plot_legend(df_positive.copy(), name, point_selector)
    chart = (
        alt.hconcat(
            scatter, legend, padding={"left": 0, "top": 0, "right": 400, "bottom": 0}
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .resolve_scale(color="shared")
    )
else:
    positive_chart = plot(
        df_positive,
        effect_x,
        effect_y,
        "Positive effects",
        point_selector,
        color_scheme,
    )
    legend_positive = plot_legend(df_positive, name, point_selector)
    negative_chart = plot(
        df_negative,
        effect_x,
        effect_y,
        "Negative effects",
        point_selector,
        color_scheme,
    )
    legend_negative = plot_legend(df_negative, name, point_selector)
    chart = (
        alt.hconcat(
            positive_chart,
            legend_positive,
            negative_chart,
            legend_negative,
            padding={"left": 0, "top": 0, "right": 400, "bottom": 0},
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .resolve_scale(color="shared")
    )

chart.save(snakemake.output[0])
