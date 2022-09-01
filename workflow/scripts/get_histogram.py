from turtle import title
import altair as alt
from altair_saver import save
import pandas as pd
import numpy as np
import pysam
from scipy.stats import gaussian_kde
from scipy import stats
import sys
import json
sys.stderr = open(snakemake.log[0], "w")


samples = snakemake.params["samples"]
f = open(snakemake.params["read_length"])
read_length = json.load(f)
f.close()


allsamp_hist_full = pd.DataFrame([])
allsamp_hist = pd.DataFrame([])
for each_sample in samples:
    sample_name = each_sample.split("/")[2]
    print(each_sample)
    print(sample_name)
    #Bam sorting and Indexing
    pysam.sort(each_sample + "/pseudoalignments.bam", "-o", each_sample + "/pseudoalignments_sort.bam")
    pysam.index("%s/pseudoalignments_sort.bam" % each_sample)
    #Bam file reading
    bam_file = pysam.AlignmentFile(each_sample + '/pseudoalignments_sort.bam',"rb")
    bam_header = bam_file.header.to_dict()
    trans_length_data = pd.DataFrame(bam_header.get('SQ'))
    trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
    #Aligned text file reading
    align_bam_txt = pd.read_csv('results/QC/' + sample_name + '.aligned.txt', 
        sep="\t",names=["read_Name", "Transcript_ID", "Start"])
    #Both transcript len and start postion are merged based on same transcript ID
    merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
    #Each read postion is calcuated
    merge_data[sample_name + '_read_dist'] = merge_data['LN'] - merge_data['Start']
    #import pdb; pdb.set_trace()
    #Python will evaluate each step
    Freq, bins = np.histogram(merge_data[sample_name + '_read_dist'], bins =read_length, range=[0,20000])
    Freq_full, bins_full = np.histogram(merge_data[sample_name + '_read_dist'], bins =read_length, range=[0,max(merge_data['LN'])])
    hist = pd.DataFrame({'sample_Name': sample_name, 'Freq': Freq, 'bins': bins[:-1]})
    hist_full = pd.DataFrame({'sample_Name': sample_name, 'Freq_full': Freq_full, 'bins_full': bins_full[:-1]})
    allsamp_hist = pd.concat([allsamp_hist, hist])
    allsamp_hist_full = pd.concat([allsamp_hist_full, hist_full])

#Histogram for first 20,000 transcript len
hist = alt.Chart(allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins', 
   scale=alt.Scale(domain=[0, 20000]), title="difference between transcript length and read start"), 
        y =alt.Y('Freq:Q', title = 'Count of Records'), 
            color = 'sample_Name').transform_filter(alt.FieldRangePredicate(field='bins', range=[0, 20000]))

allsamp_hist['read_length'] = read_length
chart_read_length = alt.Chart(allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
final_chart = (hist+chart_read_length)

#Histogram for full len transcript
hist_full = alt.Chart(allsamp_hist_full).mark_line(interpolate='step-after').encode(x = alt.X('bins_full', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_full:Q', title = 'Count of Records'), 
        color = 'sample_Name')

allsamp_hist_full['read_length'] = read_length
chart_read_length_full = alt.Chart(allsamp_hist_full).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')

final_chart_full = (hist_full+chart_read_length_full)

final_hist =alt.hconcat(final_chart, final_chart_full)

final_hist.save(snakemake.output[0])
