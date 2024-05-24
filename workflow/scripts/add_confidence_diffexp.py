import pandas as pd


def custom_sort(row, first_b_val):
    print("Row")
    print(row.head())
    print("Row")
    if row[first_b_val] < 0:
        return row[f"{first_b_val}_lower_se"]
    else:
        return row[f"{first_b_val}_upper_se"]


file_path = snakemake.input["genes_representative"]
df = pd.read_csv(file_path, sep='\t')

matching_columns = [col for col in df.columns if col.startswith(
    'b_') and not col.endswith('_se')]

for matching_column in matching_columns:

    df[f"{matching_column}_lower_se"] = df[matching_column] - \
        df[f"{matching_column}_se"]
    df[f"{matching_column}_upper_se"] = df[matching_column] + \
        df[f"{matching_column}_se"]

    df.drop(columns=[f"{matching_column}_se"], inplace=True)


print("Df")
print(df.head())
print("RowDf")


# Sort columns
b_column_order = []
for prefix in matching_columns:  # Liste der Präfixe in der gewünschten Reihenfolge
    # Liste der Suffixe in der gewünschten Reihenfolge
    for suffix in ['_lower_se', '', '_upper_se']:
        b_column_order.extend([f"{prefix}{suffix}"])
other_columns = [col for col in df.columns if not col.startswith('b_')]
df = df[other_columns + b_column_order]

# Sort rows
first_b_val = matching_columns[0]
df['sort_value'] = df.apply(lambda row: abs(
    row[f"{first_b_val}_lower_se"]) if row[first_b_val] < 0 else abs(row[f"{first_b_val}_upper_se"]), axis=1)
df = df.sort_values(by='sort_value')
df.drop(columns=["sort_value"], inplace=True)

output_file_path = snakemake.output[0]
df.to_csv(output_file_path, sep='\t', index=False)
