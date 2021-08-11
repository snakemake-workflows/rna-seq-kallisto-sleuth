from string import Template
import shutil


def main(snakemake):

    # read vega js file with template vars
    # `$url`, `$sig_level`, `$beta_column` and `$beta_se_column`
    with open(snakemake.input.spec, "rt") as f:
        spec = Template(f.read())

    tsv_path = f"{snakemake.wildcards.model}.tsv"
    shutil.copyfile(snakemake.input.tsv, snakemake.output.tsv)
    sig_level = snakemake.params.sig_level_volcano
    primary_var = snakemake.params.primary_variable

    # find column that matches primary variable
    with open(snakemake.input.tsv, "rt") as f:
        header = f.readline().strip()
        columns = header.split("\t")
        primary_cols = [
            c
            for c in columns
            if c.startswith(f"b_{primary_var}") and not c.endswith("_se")
        ]
        assert len(primary_cols) == 1
        beta_col = primary_cols[0]

    # update the spec with concrete values
    json = spec.safe_substitute(
        url=tsv_path,
        sig_level=sig_level,
        beta_column=beta_col,
        beta_se_column=beta_col + "_se",
    )

    with open(snakemake.output.json, "wt") as f:
        f.write(json)


if __name__ == "__main__":
    main(snakemake)
