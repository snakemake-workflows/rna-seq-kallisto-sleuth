import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

import warnings

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
    signed_pi_columns = [col for col in df.columns if col.startswith("signed_pi_value")]
    other_columns = [
        col
        for col in df.columns
        if not col.startswith("b_") and not col.startswith("signed_pi_value")
    ]
    return df[other_columns + b_column_order + signed_pi_columns]


def sort_rows(df):
    """Sort DataFrame by the absolute value of signed_p_value of primary variable in ascending order."""
    signed_pi_start = f"signed_pi_value_{snakemake.params['model']['primary_variable']}"
    columns_with_prefix = [col for col in df.columns if col.startswith(signed_pi_start)]

    if len(columns_with_prefix) != 1:
        warnings.warn(
            f"We can only sort by one signed_pi value column, but found {len(columns_with_prefix)}\n"
            f"respective columns with prefix '{signed_pi_start}': {columns_with_prefix}\n"
            "This usually occurs, when you have more than two levels in your condition of\n"
            "interest column.\n"
            f"We will sort by the first column found: {columns_with_prefix[0]}"
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
