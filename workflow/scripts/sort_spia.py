import pandas as pd

spia_table = snakemake.input["spia_table"]


df = pd.read_csv(spia_table, sep='\t')
df['gene_ratio'] = "(" + df['number of DE genes per pathway'].astype(
    str) + ', ' + df['number of genes on the pathway'].astype(str) + ")"
df_sorted = df.sort_values(by='Combined FDR')

output_file_path = snakemake.output[0]
df_sorted.to_csv(output_file_path, sep='\t', index=False)
