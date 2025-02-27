import pandas as pd
import altair as alt
import sys

sys.stderr = open(snakemake.log[0], "w")


def plot(df, effect_x, effect_y):
    point_selector = alt.selection_point(fields=["sample"], empty=False)

    alt.data_transformers.disable_max_rows()
    points = (
        alt.Chart(df)
        .mark_circle(size=60)
        .encode(
            alt.X(
                effect_x,
                title=effect_x,
                scale=alt.Scale(),
                axis=alt.Axis(grid=True),
            ),
            alt.Y(
                effect_y,
                title=effect_y,
                scale=alt.Scale(),
                axis=alt.Axis(grid=True),
            ),
            alt.Color(snakemake.params["color_by"], scale=alt.Scale(scheme="set1")),
            opacity=alt.value(0.5),
            tooltip=[f"sample:N", f"{effect_x}:Q", f"{effect_y}:Q"],
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
            text=alt.condition(point_selector, "sample", alt.value("")),
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
            text=alt.condition(point_selector, "sample", alt.value("")),
        )
    )

    chart = (
        alt.layer(points, text_background, text)
        .add_params(point_selector)
        .interactive()
    )
    return chart


df = pd.read_csv(snakemake.input[0], sep="\t")

pc1_pc2 = plot(df, "PC1", "PC2")
pc3_pc4 = plot(df, "PC3", "PC4")

final_chart = alt.hconcat(pc1_pc2, pc3_pc4).resolve_scale(color="shared")
final_chart.save(snakemake.output["pca"], inline=True)
