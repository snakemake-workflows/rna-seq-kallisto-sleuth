from pickle import FALSE
import altair as alt
from altair_saver import save
import os
import pandas as pd
import numpy as np
import pysam
from altair_transform import transform_chart
from pandas import DataFrame
"""
pysam.sort("-o", 'test_file/pseudoalignments_sort.bam', 'results/kallisto/SRR8309094_01C_1-1/pseudoalignments.bam')
pysam.index('test_file/pseudoalignments_sort.bam')
pure_bam = pysam.AlignmentFile('test_file/pseudoalignments_sort.bam',"rb")
for read in pure_bam:
    data = [{'read_name' : read.query_name, 'Transcript_ID' : read.reference_name, 'Start' : read.pos}]
    df = pd.DataFrame(data)
print(df)
"""
pysam.sort("-o", '../../results/kallisto/SRR8309094_01C_1-1/pseudoalignments_sort.bam', '../../results/kallisto/SRR8309094_01C_1-1/pseudoalignments.bam')
pysam.index('../../results/kallisto/SRR8309094_01C_1-1/pseudoalignments_sort.bam')
bam_file = pysam.AlignmentFile('../../results/kallisto/SRR8309094_01C_1-1/pseudoalignments_sort.bam',"rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get('SQ'))
trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
#bam_it = bam_file.fetch(until_eof = True)
#qnames = bam_file.fetch.query_name()
#qnames = [read.query_name for read in bam_it]
#qnames = np.asarray(qnames)
#print(bam_df)
#align_bam_txt = pd.read_csv('../../QC/aligned_POS_100.txt', sep="\t",names=["read_Name", "Transcript_ID", "Start"])
align_bam_txt = pd.read_csv('/home/manuel/rna-seq-kallisto-sleuth/results/kallisto/SRR8309094_01C_1-1/SRR8309094_01C_1-1_aligned.txt', sep="\t",names=["read_Name", "Transcript_ID", "Start"])
#align_bam_txt['max_rep_pos']=align_bam_txt.groupby('Transcript_ID')['Start'].transform(max)
final_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID') 
final_data['read_fract_num'] = final_data['Start']/final_data['LN']
print(final_data)
#clean_data = align_bam_txt[~align_bam_txt.isin([np.nan, np.inf, -np.inf]).any(1)]
# clean_data.to_csv('output_File2.csv')
hist = alt.Chart(final_data).mark_bar().encode(x = alt.X('read_fract_num', title="read position (relative to transcript length)",
                                                      bin = alt.BinParams(maxbins = 60)),
                                             y = 'count()')
#https://github.com/altair-viz/altair-transform 
chart = transform_chart(hist)
#chart.save(snakemake.output["histogram"])
chart.save('histogram.html')
