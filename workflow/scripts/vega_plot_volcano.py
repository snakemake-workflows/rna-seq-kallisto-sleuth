from string import Template
import pandas as pd
from io import StringIO


HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title></title>
</head>
<body>
    <div id="vis" style="width: 100%; height: 90vh;"></div>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <script>
        var vega_specs = $json;
        vegaEmbed('#vis', vega_specs);
    </script>
</body>
</html>
"""


def main(snakemake):

    # read vega js file with template vars
    # `$data`, `$sig_level`, `$beta_column` and `$beta_se_column`
    with open(snakemake.input.spec, "rt") as f:
        spec = Template(f.read())

    sig_level = snakemake.params.sig_level_volcano
    primary_var = snakemake.params.primary_variable

    # find column that matches primary variable
    df: pd.DataFrame = pd.read_csv(snakemake.input.tsv, sep="\t")

    primary_cols = [
        c
        for c in list(df.columns)
        if c.startswith(f"b_{primary_var}") and not c.endswith("_se")
    ]
    if len(primary_cols) > 1:
        print("WARNING: found {len(primary_cols)} possible primary variables")
    beta_col = primary_cols[0]

    # only keep columns needed for plot
    df = df[["ens_gene", "ext_gene", "target_id", "qval", beta_col, beta_col + "_se"]]

    # nan / NA / None values do not get plotted, so remove respective rows
    df = df.dropna()

    data = StringIO()
    df.to_csv(data, sep="\t", index=False)

    # update the spec with concrete values
    json = spec.safe_substitute(
        data=data.getvalue().replace("\t", r"\t").replace("\n", r"\n"),
        sig_level=sig_level,
        beta_column=beta_col,
        beta_se_column=beta_col + "_se",
    )

    with open(snakemake.output.json, "wt") as f:
        f.write(json)



if __name__ == "__main__":
    main(snakemake)
