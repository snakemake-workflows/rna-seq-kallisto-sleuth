import pandas as pd


def process_columns(df):
    """compute confidence interval for every column starting with b_"""
    matching_columns = [col for col in df.columns if col.startswith(
        'b_') and not col.endswith('_se')]

    for col in matching_columns:
        df[f"{col}_lower"] = df[col] - df[f"{col}_se"]
        df[f"{col}_upper"] = df[col] + df[f"{col}_se"]

    return df, matching_columns


def sort_columns(df, matching_columns):
    b_column_order = [f"{prefix}{suffix}" for prefix in matching_columns for suffix in [
        '_lower', '', '_upper', '_se']]
    other_columns = [col for col in df.columns if not col.startswith('b_')]
    return df[other_columns + b_column_order]


def sort_rows(df, primary_variable):
    """Sort DataFrame by the absolute value of signed_p_value of primary variable in ascending order."""
    df = df.reindex(
        df['signed_pi_value_' + primary_variable + '+'].abs().sort_values(descending=True).index)
    return df


# def sort_rows(df, first_b_val):
#     """Sort by b_vals if b_val < 0 sort by lower interval limit else by upper limit"""
#     df['sort_value'] = df.apply(lambda row: abs(
#         row[f"{first_b_val}_lower"]) if row[first_b_val] < 0 else abs(row[f"{first_b_val}_upper"]), axis=1)
#     df = df.sort_values(by='sort_value')
#     df.drop(columns=["sort_value"], inplace=True)
#     return df


df = pd.read_csv(snakemake.input[0], sep='\t')
df, matching_columns = process_columns(df)
df = sort_columns(df, matching_columns)
# df = sort_rows(df, matching_columns[0])
df = sort_rows(df, snakemake.params['model']['primary_variable'])
df = df.dropna(subset=matching_columns, how='all')

df.to_csv(snakemake.output[0], sep='\t', index=False)
