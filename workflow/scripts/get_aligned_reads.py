from turtle import title
import pandas as pd
import numpy as np
import pysam
import sys
#pd.set_option("display.max_rows", None, "display.max_columns", None)
samples = snakemake.params["samples"]
#samples = "results/kallisto_cds/Mel-86c_p53_Ctrl-1"
canonical_ids = pd.read_csv(snakemake.input['canonical_ids'], header=None, names=['Transcript'])
#canonical_ids = pd.read_csv("/home/manuel/paschen-quantseq-melanoma/resources/canonical_ids.csv", header=None, names=['Transcript'])
sample_name = samples.split("/")[2]
print(sample_name)
#Bam file reading
#pysam.sort(samples + "/pseudoalignments.bam", "-o", samples + "/pseudoalignments_sort.bam")
#pysam.index("%s/pseudoalignments_sort.bam" % samples)
bam_file = pysam.AlignmentFile(samples + '/pseudoalignments_sort.bam',"rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get('SQ'))
trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
print("got the transcript lengths")
#Aligned text file reading

align_bam_txt = pd.read_csv(snakemake.input['aligned_reads'], 
    sep="\t",names=["read_name","Transcript_ID", "Start", "read","Quality"])
#align_bam_txt = pd.read_csv("results/QC/Mel-86c_p53_Ctrl-1.aligned.txt", 
#    sep="\t",names=["read_name","Transcript_ID", "Start", "read","Quality"])
align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
align_bam_txt["Transcript"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[0]
print("read the aligned file")
merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
forward_strand = merge_data.loc[merge_data['Strand'] == '1']
forward_strand[sample_name + '_forward_strand'] = forward_strand['LN'] - forward_strand['Start']
aligned_reads = forward_strand.loc[forward_strand.groupby(['read_name','read'])[sample_name + '_forward_strand'].idxmin()]
reverse_strand = merge_data.loc[merge_data['Strand'] == '-1']
print("entered rev strand")
read_min = reverse_strand.loc[reverse_strand.groupby(['read_name','read'])['Start'].idxmin()]
#print(read_min)
aligned_reads = pd.concat([aligned_reads, read_min])
print("entered rev completed")
final_aligned_reads = aligned_reads.merge(canonical_ids, on="Transcript")
#final_aligned_reads.to_csv("/home/manuel/paschen-quantseq-melanoma/20220928/aligned_file.txt", sep='\t', index=False)
final_aligned_reads['read_name'].to_csv(snakemake.output['aligned_read_names'], sep='\t', index=False)