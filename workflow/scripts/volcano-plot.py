import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure
import os


def volcano(
    data: pd.DataFrame, x: str, y: str = "pval", qval: float = 0.01, title: str = ""
) -> Figure:
    assert x.startswith("b_")
    x_name = x[2:]

    data["text"] = ""
    data["significant"] = data["qval"] < qval

    fig = px.scatter(
        data,
        x=x,
        y=y,
        # text="text",
        hover_name="ext_gene",
        hover_data={
            "target_id": True,
            "pval": ":.6f",
            "qval": ":.6f",
            f"signed_pi_value_{x_name}": ":.6f",
            "text": False,
            "significant": False,
        },
        size=f"{x}_se",
        # errors_x=f"{x}_se",
        color="significant",
        color_discrete_map={True: "#34AB56", False: "#888888"},
        opacity=0.66,
        title=title,
        log_y=True,
        labels={x: "beta value", y: "q-value"},
    )
    fig.add_hline(
        y=qval,
        line_dash="dash",
        annotation_text=f"qval: {qval}",
    )
    fig.update_layout(
        template="plotly_white", showlegend=True, yaxis=dict(autorange="reversed")
    )
    return fig


def read_diffexp_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # to dropna or to fillna, that is the question
    df.dropna(inplace=True)
    return df


def main(snakemake):
    model = snakemake.wildcards.model
    os.makedirs(snakemake.output.volcano_plot)
    data = read_diffexp_matrix(snakemake.input.tsv)
    covariates = [
        c[2:] for c in data.columns if c.startswith("b_") and not c.endswith("_se")
    ]
    for covariate in covariates:
        sig_level = snakemake.params.sig_level_volcano
        fig = volcano(
            data,
            x="b_" + covariate,
            y="qval",
            qval=sig_level,
            title=f"volcano plot for {model} {covariate} (qval={sig_level})",
        )
        fig.write_html(os.path.join(snakemake.output.volcano_plot, covariate + ".html"))


if __name__ == "__main__":
    main(snakemake)
