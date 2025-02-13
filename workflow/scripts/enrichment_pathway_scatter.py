import pandas as pd
import altair as alt
import sys

sys.stderr = open(snakemake.log[0], "w")


def plot(df_plot, effect_x, effect_y, title, selector, color_scheme):
    alt.data_transformers.disable_max_rows()
    points = (
        alt.Chart(df_plot, title=title)
        .transform_calculate(
            order_color=f"indexof({selector.name}.Name || [], datum.Name)",
            order_shape=f"indexof({selector.name}.Name || [], datum.Name) % 5",
        )
        .mark_point(size=40)
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
                    "order_color:N", scale=alt.Scale(range=color_scheme), legend=None
                ),
                alt.value("lightgray"),
            ),
            shape=alt.condition(
                selector,
                alt.Shape(
                    "order_shape:N",
                    scale=alt.Scale(
                        domain=[0, 1, 2, 3, 4],
                        range=[
                            "triangle-up",
                            "triangle-down",
                            "square",
                            "diamond",
                            "cross",
                        ],
                    ),
                    legend=None,
                ),
                alt.value("circle"),
            ),
            tooltip=[f"{name}:N", f"{effect_x}:Q", f"{effect_y}:Q"],
        )
        .add_params(selector)
    )
    return points


def plot_legend(df_plot, name, selector, color_scheme):
    legend = (
        alt.Chart(df_plot)
        .transform_calculate(
            order_color=f"indexof({selector.name}.Name || [], datum.Name)",
            order_shape=f"indexof({selector.name}.Name || [], datum.Name) % 5",
        )
        .mark_point(size=60)
        .encode(
            alt.Y(
                name,
                axis=alt.Axis(
                    title="", ticks=False, domain=False, labelLimit=700, orient="right"
                ),
            ),
            color=alt.condition(
                selector,
                alt.Color(
                    "order_color:N", scale=alt.Scale(range=color_scheme), legend=None
                ),
                alt.value("lightgrey"),
            ),
            shape=alt.condition(
                selector,
                alt.Shape(
                    "order_shape:N",
                    scale=alt.Scale(
                        domain=[0, 1, 2, 3, 4],
                        range=[
                            "triangle-up",
                            "triangle-down",
                            "square",
                            "diamond",
                            "cross",
                        ],
                    ),
                    legend=None,
                ),
                alt.value("circle"),
            ),
        )
        .transform_filter(selector)
    )
    return legend


color_scheme = [
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]

df = pd.read_csv(snakemake.input[0], sep="\t")
df = df.head(50)
name = snakemake.params["name"]
effect_x = snakemake.params["effect_x"]
effect_y = snakemake.params["effect_y"]

df_positive = df[df[f"{effect_x}"] >= 0]
df_negative = df[df[f"{effect_x}"] < 0]
df_negative[f"{effect_x}"] = df_negative[f"{effect_x}"].abs()

point_selector = alt.selection_point(fields=["Name"], empty=False)

if df_negative.empty:
    scatter = plot(
        df_positive, effect_x, effect_y, "", point_selector, color_scheme, name
    )
    # Important: You need to copy the df in order to have different datasets, else vega does not bind the plot to a dataset
    legend = plot_legend(df_positive.copy(), name, point_selector)
    chart = (
        alt.hconcat(
            scatter, legend, padding={"left": 0, "top": 0, "right": 400, "bottom": 0}
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .resolve_scale(color="shared", shape="shared")
    )
else:
    positive_scatter = plot(
        df_positive,
        effect_x,
        effect_y,
        "Positive effects",
        point_selector,
        color_scheme,
    )
    legend_positive = plot_legend(df_positive, name, point_selector, color_scheme)
    positive_chart = alt.hconcat(
        positive_scatter,
        legend_positive,
    )
    negative_scatter = plot(
        df_negative,
        effect_x,
        effect_y,
        "Negative effects",
        point_selector,
        color_scheme,
    )
    legend_negative = plot_legend(df_negative, name, point_selector, color_scheme)
    negative_chart = alt.hconcat(
        negative_scatter,
        legend_negative,
    )

    chart = (
        alt.vconcat(
            positive_chart,
            negative_chart,
            padding={"left": 0, "top": 0, "right": 300, "bottom": 0},
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .resolve_scale(color="shared", shape="shared")
    )

chart.save(snakemake.output[0])
