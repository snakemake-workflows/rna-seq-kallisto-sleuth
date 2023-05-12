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
transcript_ids = snakemake.params["each_transcript"]

samples = snakemake.params["samples"]

f = open(snakemake.params["read_length"])
read_length = json.load(f)
f.close()
read_length = 43
fwrd_allsamp_hist = pd.DataFrame([])
fwrd_allsamp_hist_trim = pd.DataFrame([])
rev_allsamp_hist = pd.DataFrame([])
rev_allsamp_hist_trim = pd.DataFrame([])
fwrd_flag = 0
rev_flag = 0
empty_histogram =pd.DataFrame([])


# taking each sample aligned position
for each_sample in samples:
    sample_name = each_sample.split("/")[2]
    # Bam file reading
    bam_file = pysam.AlignmentFile('results/ind_transcripts/' + sample_name + '-pseudoalignments.sorted.bam',"rb")
    bam_header = bam_file.header.to_dict()
    get_trans_length = pd.DataFrame(bam_header.get('SQ'))
    get_trans_length.rename(columns={'SN': 'Transcript_ID'}, inplace=True)

    # Aligned text file reading
    align_bam_txt = pd.read_csv('results/QC/' + sample_name + '.aligned.txt', 
    sep="\t",names=["read_Name","Transcript_ID", "Start","reads", "Quality"])
    align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
    align_bam_txt["Transcript"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[0]

    # merging aligned bam text file based on transcript id from bam file
    merge_align_trans_data = align_bam_txt.merge(get_trans_length, on='Transcript_ID')
    filtered_transcript_data = merge_align_trans_data.query("Transcript == @transcript_ids")

    # Forward strand
    if ((filtered_transcript_data['Strand']=='1')).any():
        fwrd_flag = 1          
        filtered_transcript_data[sample_name + '_forward_strand'] = filtered_transcript_data['LN'] - filtered_transcript_data['Start']
        aligned_reads = filtered_transcript_data.loc[filtered_transcript_data.groupby(['read_Name','reads'])[sample_name + '_forward_strand'].idxmin()]
        Freq_fwrd, bins_fwrd = np.histogram(aligned_reads[sample_name + '_forward_strand'], bins = read_length, range=[0,max(aligned_reads['LN'])])
        hist_fwrd = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd, 'bins_foward': bins_fwrd[:-1]})
        fwrd_allsamp_hist = pd.concat([fwrd_allsamp_hist, hist_fwrd])
    elif ((filtered_transcript_data['Strand']=='-1')).any():
    # reverse strand
        rev_flag = 1
        reverse_strand = filtered_transcript_data.loc[filtered_transcript_data['Strand'] == '-1']
        read_min = reverse_strand.loc[reverse_strand.groupby(['read_Name','reads'])['Start'].idxmin()]
        Freq_rev, bins_rev = np.histogram(read_min['Start'], bins =read_length, range=[0, max(read_min['LN'])])
        hist_rev = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev,'bins_rev': bins_rev[:-1]})
        rev_allsamp_hist = pd.concat([rev_allsamp_hist, hist_rev])

if fwrd_flag == 1:
    # Histogram for full len transcript
    hist_fwrd_full = alt.Chart(fwrd_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_foward', 
    title = "difference between transcript length and read start"), y = alt.Y('Freq_forward:Q', title = 'Count of Records'), 
    color = 'sample_Name').properties(title = "forward strand transcripts (full length)")
    fwrd_allsamp_hist['read_length'] = read_length
    cht_rd_len_fwd_full = alt.Chart(fwrd_allsamp_hist).mark_rule(color = 'black',strokeDash = [3,5]).encode(
    x = 'read_length')
    Fwd_chart = (hist_fwrd_full+cht_rd_len_fwd_full)
    Fwd_chart_full = alt.hconcat(Fwd_chart)
    
    Fwd_chart_full.save(snakemake.output[0])

elif rev_flag == 1:
    # Histogram for reverse strand
    # Histogram for full len transcript
    hist_rev_full = alt.Chart(rev_allsamp_hist).mark_line(interpolate='step-after').encode(x = alt.X('bins_rev', 
    title="read start position in the transcript"), y =alt.Y('Freq_rev:Q', title = 'Count of Records'), 
    color = 'sample_Name').properties(title="Reverse strand transcripts (full length)")
    rev_allsamp_hist['read_length'] = read_length
    cht_rd_len_rev = alt.Chart(rev_allsamp_hist).mark_rule(color='black',strokeDash=[3,5]).encode(
    x='read_length')
    rev_chart =(hist_rev_full+cht_rd_len_rev)
    Rev_chart_full =alt.hconcat(rev_chart)
    
    Rev_chart_full.save(snakemake.output[0])
else:
    
    # Configure empty histogram
    empty_histogram = alt.Chart().mark_text(text='no reads aligned',size=20)
    empty_histogram.save(snakemake.output[0])