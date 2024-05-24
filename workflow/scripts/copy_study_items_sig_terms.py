import pandas as pd

file_path_enr = snakemake.input["enrichment"]
file_path_sig = snakemake.input["significant_terms"]

df_enr = pd.read_csv(file_path_enr, sep='\t')
df_sig = pd.read_csv(file_path_sig, sep='\t')

# Some entries are only in one of both tables
common_ids = df_sig[df_sig['GO'].isin(df_enr['GO'])]['GO']

# Filter on common ids
df_enr_filtered = df_enr[df_enr['GO'].isin(common_ids)]
df_sig_filtered = df_sig[df_sig['GO'].isin(common_ids)]

df_enr_filtered['study_items_sig_terms'] = df_enr_filtered['GO'].map(
    df_sig_filtered.set_index('GO')['study_items'])

# Define function to sort each group by p_uncorrected ascending


def sort_group(group):
    return group.sort_values(by='p_uncorrected', ascending=True)


# Group by 'class' and apply sorting function
df_enr_filtered_sorted = df_enr_filtered.groupby(
    'class', group_keys=False).apply(sort_group)

output_file_path = snakemake.output[0]
df_enr_filtered_sorted.to_csv(output_file_path, sep='\t', index=False)
