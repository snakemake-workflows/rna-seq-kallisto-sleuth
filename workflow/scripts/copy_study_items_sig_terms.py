import pandas as pd

file_path_enr = snakemake.input["enrichment"]
file_path_sig = snakemake.input["significant_terms"]

df_enr = pd.read_csv(file_path_enr, sep='\t')
df_sig = pd.read_csv(file_path_sig, sep='\t')

df_enr['study_items_sig_terms'] = df_sig['study_items']

output_file_path = snakemake.output[0]
df_enr.to_csv(output_file_path, sep='\t', index=False)
