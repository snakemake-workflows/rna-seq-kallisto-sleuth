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

fwrd_allsamp_hist = pd.DataFrame([])
fwrd_allsamp_hist_trim = pd.DataFrame([])
rev_allsamp_hist = pd.DataFrame([])
rev_allsamp_hist_trim = pd.DataFrame([])
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
        sep="\t",names=["read_Name","Transcript_ID", "Start","reads", "Quality"])
    align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
    #Both transcript len and start postion are merged based on same transcript ID
    merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
    #Each read postion is calcuated
    forward_strand = merge_data.loc[merge_data['Strand'] == '1']
    forward_strand[sample_name + '_forward_strand'] = forward_strand['LN'] - forward_strand['Start']
    aligned_reads = forward_strand.loc[forward_strand.groupby(['read_Name','reads'])[sample_name + '_forward_strand'].idxmin()]

    Freq_fwrd, bins_fwrd = np.histogram(aligned_reads[sample_name + '_forward_strand'], bins =read_length, range=[0,max(aligned_reads['LN'])])
    Freq_fwrd_trim, bins_fwrd_trim = np.histogram(forward_strand[sample_name + '_forward_strand'], bins =read_length, range=[0, 20000])

    hist_fwrd = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd, 'bins_foward': bins_fwrd[:-1]})
    hist_fwrd_trim = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd_trim, 'bins_foward': bins_fwrd_trim[:-1]})

    fwrd_allsamp_hist = pd.concat([fwrd_allsamp_hist, hist_fwrd])
    fwrd_allsamp_hist_trim = pd.concat([fwrd_allsamp_hist_trim, hist_fwrd_trim])
    #reverse strand
    reverse_strand = merge_data.loc[merge_data['Strand'] == '-1']
    print("entered rev strand")
    read_min = reverse_strand.loc[reverse_strand.groupby(['read_Name','reads'])['Start'].idxmin()]

    Freq_rev, bins_rev = np.histogram(read_min['Start'], bins =read_length, range=[0,max(read_min['LN'])])
    Freq_rev_trim, bins_rev_trim = np.histogram(read_min['Start'], bins =read_length, range=[0, 20000])

    hist_rev = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev, 'bins_rev': bins_rev[:-1]})
    hist_rev_trim = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev_trim, 'bins_rev': bins_rev_trim[:-1]})

    rev_allsamp_hist = pd.concat([rev_allsamp_hist, hist_rev])
    rev_allsamp_hist_trim = pd.concat([rev_allsamp_hist_trim, hist_rev_trim])
    print("entered rev completed")

#Histogram for forward strand
#Histogram for full len transcript
hist_fwrd_full = alt.Chart(fwrd_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_foward', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_forward:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the forward strand full transcript length")
#Histogram upto 20000 bp transcript length
hist_fwrd_trimd = alt.Chart(fwrd_allsamp_hist_trim).mark_line(interpolate='step-after').encode(x = alt.X('bins_foward', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_forward:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the forward strand up to 20000 base pair transcript length")
fwrd_allsamp_hist['read_length'] = read_length
fwrd_allsamp_hist_trim['read_length'] = read_length
cht_rd_len_fwd_full = alt.Chart(fwrd_allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
cht_rd_len_fwd_trim = alt.Chart(fwrd_allsamp_hist_trim).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')

#Histogram for forward strand
#Histogram for full len transcript
hist_rev_full = alt.Chart(rev_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_rev', 
    title="read start position in the transcript"), y =alt.Y('Freq_rev:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the reverse strand full transcript length")
hist_rev_trim = alt.Chart(rev_allsamp_hist_trim).mark_line(interpolate='step-after').encode(x = alt.X('bins_rev', 
    title="read start position in the transcript"), y =alt.Y('Freq_rev:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the reverse strand up to 20000 base pair transcript length")

rev_allsamp_hist['read_length'] = read_length
rev_allsamp_hist_trim['read_length'] = read_length

cht_rd_len_rev = alt.Chart(rev_allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
cht_rd_len_rev_trim = alt.Chart(rev_allsamp_hist_trim).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')


Fwd_chart_full = (hist_fwrd_full+cht_rd_len_fwd_full)
Fwd_chart_trim = (hist_fwrd_trimd+cht_rd_len_fwd_trim)

rev_chart_full =(hist_rev_full+cht_rd_len_rev)
rev_chart_trim =(hist_rev_trim+cht_rd_len_rev_trim)

Fwd_chart =alt.hconcat(Fwd_chart_trim, Fwd_chart_full)
Rev_chart =alt.hconcat(rev_chart_trim, rev_chart_full)

Final_chart = alt.vconcat(Fwd_chart, Rev_chart)
Final_chart.save(snakemake.output[0])
