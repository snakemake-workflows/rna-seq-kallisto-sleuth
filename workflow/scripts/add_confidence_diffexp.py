import pandas as pd

file_path = snakemake.input["genes_representative"]

df = pd.read_csv(file_path, sep='\t')

matching_column = [col for col in df.columns if col.startswith('b_') and not col.endswith(
    '_part2') and not col.endswith('_se')][0]

df[f"{matching_column}_lower_se"] = df[matching_column] - \
    df[f"{matching_column}_se"]
df[f"{matching_column}_upper_se"] = df[matching_column] + \
    df[f"{matching_column}_se"]

output_file_path = snakemake.output[0]
df.to_csv(output_file_path, sep='\t', index=False)
