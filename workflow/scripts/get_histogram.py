from pickle import FALSE
import altair as alt
from altair_saver import save
import os
import pandas as pd
import numpy as np
import pysam
from altair_transform import transform_chart
from pandas import DataFrame

path = snakemake.input["bam_path"]
pysam.sort(path + "/pseudoalignments.bam", "-o", path + "/pseudoalignments_sort.bam")
pysam.index("%s/pseudoalignments_sort.bam" % path)


bam_file = pysam.AlignmentFile(path + '/pseudoalignments_sort.bam',"rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get('SQ'))
trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)

align_bam_txt = pd.read_csv(snakemake.input["aligned_file"], sep="\t",names=["read_Name", "Transcript_ID", "Start"])
final_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID') 
final_data['read_fract_num'] = final_data['Start']/final_data['LN']
print(final_data)
# clean_data.to_csv('output_File2.csv')
hist = alt.Chart(final_data).mark_bar().encode(x = alt.X('read_fract_num', title="read position (relative to transcript length)",
                                                      bin = alt.BinParams(maxbins = 60)),
                                             y = 'count()')
#https://github.com/altair-viz/altair-transform 
chart = transform_chart(hist)
chart.save(snakemake.output["histogram"])
