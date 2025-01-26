import pandas as pd
import altair as alt

df = pd.read_csv(snakemake.input[0], sep="\t")

name = snakemake.params["name"]
title = snakemake.params["title"]
effect_x = snakemake.params["effect_x"]
effect_y = snakemake.params["effect_y"]
df_top_50 = df.head(50)

point_selector = alt.selection_single(fields=[f"{name}"], empty=False)


# Filter out rows where either effect_x or effect_y is zero because of logarithmic scale
# point_selector = alt.selection_single(fields=[f"{name}"], empty=False)

alt.data_transformers.disable_max_rows()
points = (
    alt.Chart(df_top_50)
    .mark_circle(
        size=50,
        # tooltip=[f"{name}:N", f"{effect_x}:Q", f"{effect_y}:Q"],
    )
    .encode(
        alt.X(
            effect_x,
            title=effect_x,
            scale=alt.Scale(type="log", nice=False),
            axis=alt.Axis(grid=True),
        ),
        alt.Y(
            effect_y,
            title=effect_y,
            scale=alt.Scale(type="log", nice=False),
            axis=alt.Axis(grid=True),
        ),
        tooltip=[f"{name}:N", f"{effect_x}:Q", f"{effect_y}:Q"],
        color=alt.condition(
            point_selector,
            alt.value("red"),  # Selected points are red
            alt.value("lightgray"),  # Default points are light gray
        ),
        #         opacity=alt.condition(
        #             point_selector, alt.value(1), alt.value(0.2)
        #         ),
    )
)


text_background = (
    alt.Chart(df_top_50)
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
        text=alt.condition(point_selector, f"{name}", alt.value("")),
    )
)

text = (
    alt.Chart(df_top_50)
    .mark_text(
        align="left",
        baseline="middle",
        dx=5,
        dy=-5,
    )
    .encode(
        x=effect_x,
        y=effect_y,
        text=alt.condition(point_selector, f"{name}", alt.value("")),
    )
)

chart = (
    # alt.layer(line, points, text_background, text)
    alt.layer(points, text_background, text)
    .add_params(point_selector)
    .properties(title=title)
    .interactive()
)


# main_chart = (
#     alt.Chart(df_top_50)
#     .mark_circle(size=60)
#     .encode(
#         alt.Y(
#             effect_y,
#             title=effect_y,
#             scale=alt.Scale(type="log"),
#             axis=alt.Axis(grid=True),
#         ),
#         alt.X(
#             effect_x,
#             title=effect_x,
#             scale=alt.Scale(type="log"),
#             axis=alt.Axis(grid=True),
#         ),
#         color=f"{name}",
#         tooltip=[f"{name}:N", f"{effect_x}:Q", f"{effect_y}:Q"],
#         opacity=alt.condition(
#             point_selector, alt.value(1), alt.value(0.2)
#         ),  # Opacity ändern je nach Auswahl
#     )
#     .add_params(point_selector)  # Parameter hinzufügen, die die Auswahl steuern
# )

chart.save(snakemake.output[0], inline=True)
