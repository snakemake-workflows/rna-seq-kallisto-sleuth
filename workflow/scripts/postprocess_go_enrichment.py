import pandas as pd
from ast import literal_eval as make_tuple


def calculate_enrichment(ratio_in_study_str, ratio_in_pop_str):
    """Calculate enrichment based on ratios provided as strings."""
    ratio_in_study = make_tuple(ratio_in_study_str)
    ratio_in_pop = make_tuple(ratio_in_pop_str)
    enrichment_study = ratio_in_study[0] / ratio_in_study[1]
    enrichment_pop = ratio_in_pop[0] / ratio_in_pop[1]
    return enrichment_study, enrichment_pop


def sort_group(group):
    """Sort a DataFrame group by the 'p_uncorrected' column in ascending order."""
    return group.sort_values(by='p_uncorrected', ascending=True)


# Load data
df_enr = pd.read_csv(snakemake.input["enrichment"], sep='\t')
df_sig = pd.read_csv(snakemake.input["significant_terms"], sep='\t')

# TODO Would be better with merge but does not work right now
# # Merge the two dataframes on the 'GO' column to keep only common GO terms
# df_merged = df_sig.join(df_enr.set_index('GO'), on='GO', rsuffix='_enr')
# # Add study items from significant terms to the merged dataset
# df_merged['study_items_sig_terms'] = df_merged['study_items']
# Only keep data if GO term exists in both tables
common_ids = df_sig[df_sig['GO'].isin(df_enr['GO'])]['GO']
df_enr_filtered = df_enr[df_enr['GO'].isin(common_ids)]
df_sig_filtered = df_sig[df_sig['GO'].isin(common_ids)]

# Add study items from significant terms to dataset
df_enr_filtered['study_items_sig_terms'] = df_enr_filtered['GO'].map(
    df_sig_filtered.set_index('GO')['study_items'])

# Sort and calculate enrichment ratios
df_merged_sorted = df_enr_filtered.groupby(
    'class', group_keys=False).apply(sort_group)

if not df_merged_sorted.empty:
    df_merged_sorted['enrichment'] = df_merged_sorted.apply(
        lambda row: calculate_enrichment(row['ratio_in_study'], row['ratio_in_pop']), axis=1)
else:
    df_merged_sorted['enrichment'] = None

# Save the result to a file
df_merged_sorted.to_csv(snakemake.output[0], sep='\t', index=False)
