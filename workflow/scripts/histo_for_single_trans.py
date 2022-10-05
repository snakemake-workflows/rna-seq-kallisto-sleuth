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
pd.set_option("display.max_rows", None, "display.max_columns", None)
transcript_ids = snakemake.params["transcripts"]
#samples = ["Mel-86c_p49_Ctrl-1", "Mel-86c_p51_Ctrl-1", "Mel-86c_p53_Ctrl-1", "Mel-86c_p49_IFNg50-1", "Mel-86c_p51_IFNg50-1", "Mel-86c_p53_IFNg50-1", "Mel-86c_p49_Plx-1", "Mel-86c_p51_Plx-1", "Mel-86c_p53_Plx-1", "Mel-86c_p49_IFNg50_Plx-1", "Mel-86c_p51_IFNg50_Plx-1", "Mel-86c_p53_IFNg50_Plx-1"]
#samples = ["10435_17-1","15281_17_A-1", "1697_18-1","1738_17-1","1797_18_b-1","1802193_a-1"]
#samples = ["10435_17-1"]
samples = ["Mel-86c_p49_Ctrl-1"]
#samples = ["Mel-86c_p49_IFNg50-1"]
read_length= 43
fwrd_allsamp_hist = pd.DataFrame([])
fwrd_allsamp_hist_trim = pd.DataFrame([])
rev_allsamp_hist = pd.DataFrame([])
rev_allsamp_hist_trim = pd.DataFrame([])
print(transcript_ids)
"""
for each_sample in samples:
    #sample_name = each_sample.split("/")[4]
    print(each_sample)
    #print(sample_name)
    #Bam sorting and Indexing
    #pysam.sort("/home/manuel/paschen-quantseq-melanoma/results/kallisto_cds/" + each_sample + "/pseudoalignments.bam", "-o", "/home/manuel/paschen-quantseq-melanoma/results/kallisto_cds/" + each_sample + "/pseudoalignments_sort.bam")
    #pysam.index("%s/pseudoalignments_sort.bam" % "/home/manuel/paschen-quantseq-melanoma/results/kallisto_cds" + each_sample )
    #Bam file reading
    bam_file = pysam.AlignmentFile('/home/manuel/paschen-quantseq-melanoma/results/kallisto_cds/' + each_sample + '/pseudoalignments_sort.bam',"rb")
    bam_header = bam_file.header.to_dict()
    trans_length_data = pd.DataFrame(bam_header.get('SQ'))
    trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
    print("got the transcript lengths")
    #Aligned text file reading


    align_bam_txt = pd.read_csv('/home/manuel/paschen-quantseq-melanoma/ENST00000473504.1_ENST00000388918.10/' + each_sample + '.aligned.txt', 
        sep=" ",names=["Transcript_ID", "Start", "read"])
    align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
    print("read the aligned file")
    merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
    #print(merge_data)
    #sel_transcript=df.loc[df['column_name'] == some_value]
    forward_strand = merge_data.loc[merge_data['Strand'] == '1']
    fwd_read_min = forward_strand.loc[forward_strand.groupby('read')['Start'].idxmin()]
    fwd_read_min[each_sample + '_forward_strand'] = fwd_read_min['LN'] - fwd_read_min['Start']
    #print(forward_strand)
    trans_diff_min = fwd_read_min.loc[fwd_read_min.groupby('read')[each_sample + '_forward_strand'].idxmin()]
    #print(trans_diff_min)
    Freq_fwrd, bins_fwrd = np.histogram(trans_diff_min[each_sample + '_forward_strand'], bins =43, range=[0,max(trans_diff_min['LN'])])
    #Freq_fwrd_trim, bins_fwrd_trim = np.histogram(read_max[each_sample + '_forward_strand'], bins =60, range=[0, 20000])
    hist_fwrd = pd.DataFrame({'sample_Name': each_sample, 'Freq_forward': Freq_fwrd, 'bins_foward': bins_fwrd[:-1]})
    #hist_fwrd_trim = pd.DataFrame({'sample_Name': each_sample, 'Freq_forward': Freq_fwrd_trim, 'bins_foward': bins_fwrd_trim[:-1]})
    fwrd_allsamp_hist = pd.concat([fwrd_allsamp_hist, hist_fwrd])
    #fwrd_allsamp_hist_trim = pd.concat([fwrd_allsamp_hist_trim, hist_fwrd_trim])
    
    #reverse strand
    reverse_strand = merge_data.loc[merge_data['Strand'] == '-1']
    print("entered rev strand")
    read_min = reverse_strand.loc[reverse_strand.groupby('read_name')['Start'].idxmin()]
    print(read_min)
    Freq_rev, bins_rev = np.histogram(read_min['Start'], bins =60, range=[0,max(read_min['LN'])])
    Freq_rev_trim, bins_rev_trim = np.histogram(read_min['Start'], bins =60, range=[0, 20000])
    hist_rev = pd.DataFrame({'sample_Name': each_sample, 'Freq_rev': Freq_rev, 'bins_rev': bins_rev[:-1]})
    hist_rev_trim = pd.DataFrame({'sample_Name': each_sample, 'Freq_rev': Freq_rev_trim, 'bins_rev': bins_rev_trim[:-1]})
    rev_allsamp_hist = pd.concat([rev_allsamp_hist, hist_rev])
    rev_allsamp_hist_trim = pd.concat([rev_allsamp_hist_trim, hist_rev_trim])
    print("entered rev completed")
    
print(trans_diff_min)
print(fwrd_allsamp_hist)
#Histogram for forward strand
#Histogram for full len transcript
hist_fwrd_full = alt.Chart(fwrd_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_foward', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_forward:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the forward strand full transcript length")
#Histogram upto 20000 bp transcript length
#hist_fwrd_trimd = alt.Chart(fwrd_allsamp_hist_trim).mark_line(interpolate='step-after').encode(x = alt.X('bins_foward', 
#    title="difference between transcript length and read start"), y =alt.Y('Freq_forward:Q', title = 'Count of Records'), 
#        color = 'sample_Name').properties(title="QC plot of the forward strand up to 20000 base pair transcrip length")
fwrd_allsamp_hist['read_length'] = read_length
#fwrd_allsamp_hist_trim['read_length'] = read_length
cht_rd_len_fwd_full = alt.Chart(fwrd_allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
#cht_rd_len_fwd_trim = alt.Chart(fwrd_allsamp_hist_trim).mark_rule(color='black',strokeDash=[3,5]).encode(
#   x='read_length')

#Histogram for forward strand
#Histogram for full len transcript

hist_rev_full = alt.Chart(rev_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_rev', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_rev:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the reverse strand full transcript length")
hist_rev_trim = alt.Chart(rev_allsamp_hist_trim).mark_line(interpolate='step-after').encode(x = alt.X('bins_rev', 
    title="difference between transcript length and read start"), y =alt.Y('Freq_rev:Q', title = 'Count of Records'), 
        color = 'sample_Name').properties(title="QC plot of the reverse strand up to 20000 base pair transcrip length")

rev_allsamp_hist['read_length'] = read_length
rev_allsamp_hist_trim['read_length'] = read_length

cht_rd_len_rev = alt.Chart(rev_allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
cht_rd_len_rev_trim = alt.Chart(rev_allsamp_hist_trim).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')


Fwd_chart_full = (hist_fwrd_full+cht_rd_len_fwd_full)
#Fwd_chart_trim = (hist_fwrd_trimd+cht_rd_len_fwd_trim)

#rev_chart_full =(hist_rev_full+cht_rd_len_rev)
#rev_chart_trim =(hist_rev_trim+cht_rd_len_rev_trim)

#Fwd_chart =alt.hconcat(Fwd_chart_trim, Fwd_chart_full)
Fwd_chart =alt.hconcat( Fwd_chart_full)
#Rev_chart =alt.hconcat(rev_chart_trim, rev_chart_full)


print("concatination completed")
#Final_chart = alt.vconcat(Fwd_chart, Rev_chart)
Final_chart = alt.vconcat(Fwd_chart)
Final_chart.save("/home/manuel/paschen-quantseq-melanoma/ENST00000473504.1_ENST00000388918.10/Mel-86c_p49_IFNg50-1_ENST00000388918.10_QC_plot.html")
"""