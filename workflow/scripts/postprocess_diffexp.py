import pandas as pd


def process_columns(df):
    """compute confidence interval for every column starting with b_"""
    matching_columns = [
        col for col in df.columns if col.startswith("b_") and not col.endswith("_se")
    ]

    for col in matching_columns:
        df[f"{col}_lower"] = df[col] - df[f"{col}_se"]
        df[f"{col}_upper"] = df[col] + df[f"{col}_se"]

    return df, matching_columns


def sort_columns(df, matching_columns):
    b_column_order = [
        f"{prefix}{suffix}"
        for prefix in matching_columns
        for suffix in ["_lower", "", "_upper", "_se"]
    ]
    other_columns = [col for col in df.columns if not col.startswith("b_")]
    return df[other_columns + b_column_order]


def sort_rows(df):
    """Sort DataFrame by the absolute value of signed_p_value of primary variable in ascending order."""
    signed_pi_start = f"signed_pi_value_{snakemake.params['model']['primary_variable']}"
    columns_with_prefix = [col for col in df.columns if col.startswith(signed_pi_start)]

    if len(columns_with_prefix) != 1:
        raise ValueError(
            f"Expected exactly one column starting with '{signed_pi_start}', found {len(columns_with_prefix)}"
        )

    signed_pi_col = columns_with_prefix[0]

    df_sorted = df.sort_values(by=signed_pi_col, key=lambda x: x.abs(), ascending=False)
    return df_sorted


df = pd.read_csv(snakemake.input[0], sep="\t")
df, matching_columns = process_columns(df)
df = sort_columns(df, matching_columns)
df = sort_rows(df)
df = df.dropna(subset=matching_columns, how="all")
df.to_csv(snakemake.output[0], sep="\t", index=False)
