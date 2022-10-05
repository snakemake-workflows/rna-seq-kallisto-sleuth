import pandas as pd
import numpy as np
import pysam
import sys

sys.stderr = open(snakemake.log[0], "w")

samples = snakemake.params["samples"]
canonical_ids = pd.read_csv(snakemake.input['canonical_ids'], header=None, names=['Transcript'])
sample_name = samples.split("/")[2]
#Bam file reading
pysam.sort(samples + "/pseudoalignments.bam", "-o", samples + "/pseudoalignments_sort.bam")
pysam.index("%s/pseudoalignments_sort.bam" % samples)
bam_file = pysam.AlignmentFile(samples + '/pseudoalignments_sort.bam',"rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get('SQ'))
trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
#Aligned text file reading
align_bam_txt = pd.read_csv(snakemake.input['aligned_reads'], 
    sep="\t",names=["read_name","Transcript_ID", "Start", "read","Quality"])
align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
align_bam_txt["Transcript"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[0]
merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
#reads aligned to forward strand
forward_strand = merge_data.loc[merge_data['Strand'] == '1']
forward_strand[sample_name + '_forward_strand'] = forward_strand['LN'] - forward_strand['Start']
aligned_reads = forward_strand.loc[forward_strand.groupby(['read_name','read'])[sample_name + '_forward_strand'].idxmin()]
#reads aligned to reverse strand
reverse_strand = merge_data.loc[merge_data['Strand'] == '-1']
read_min = reverse_strand.loc[reverse_strand.groupby(['read_name','read'])['Start'].idxmin()]
aligned_reads = pd.concat([aligned_reads, read_min])
final_aligned_reads = aligned_reads.merge(canonical_ids, on="Transcript")
final_aligned_reads['read_name'].to_csv(snakemake.output['aligned_read_names'], sep='\t', index=False)