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
import csv 

sys.stderr = open(snakemake.log[0], "w")

samples = snakemake.input["aligned_file"]
f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()

transcript_ids = snakemake.params["each_transcript"]

fwrd_allsamp_hist = pd.DataFrame([])
fwrd_allsamp_hist_trim = pd.DataFrame([])
fwrd_allsamp_hist_fil = pd.DataFrame([])
rev_allsamp_hist = pd.DataFrame([])
rev_allsamp_hist_trim = pd.DataFrame([])
rev_allsamp_hist_fil = pd.DataFrame([])
sample_file = samples.split("/")[2]
sample_name = sample_file.split(".aligned.txt")[0]

#Bam file reading
bam_file = pysam.AlignmentFile(snakemake.input["samtools_sort"],"rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get('SQ'))
trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)

#Aligned text file reading
align_bam_txt = pd.read_csv(snakemake.input["aligned_file"], 
    sep="\t",names=["read_Name","Transcript_ID", "Start","reads", "Quality"])
align_bam_txt["Strand"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[1]
align_bam_txt["Transcript"] =align_bam_txt['Transcript_ID'].str.split('_', 1).str[0]

#Both transcript len and start postion are merged based on same transcript ID
merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
#filtered_transcript_data = merge_data.query("Transcript == @transcript_ids")

#Forward strand
forward_strand = merge_data.loc[merge_data['Strand'] == '1']

#Each read postion is calcuated
forward_strand[sample_name + '_forward_strand'] = forward_strand['LN'] - forward_strand['Start']
aligned_reads = forward_strand.loc[forward_strand.groupby(['read_Name','reads'])[sample_name + '_forward_strand'].idxmin()]

if aligned_reads['Transcript'].str.contains(transcript_ids).any():
    print("entered single transcript")
    fwrd_filtered_transcript_data = aligned_reads.query("Transcript == @transcript_ids")
    Freq_fwrd_fil, bins_fwrd_fil = np.histogram(aligned_reads[sample_name + '_forward_strand'], bins =read_length, range=[0,max(aligned_reads['LN'])])
    hist_fwrd_fil = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd_fil, 'bins_foward': bins_fwrd_fil[:-1]})
    fwrd_allsamp_hist_fil = pd.concat([fwrd_allsamp_hist_fil, hist_fwrd_fil])
elif transcript_ids == "all":
    print("entered all transcript")
    #Values are added to corresponding bins
    Freq_fwrd, bins_fwrd = np.histogram(aligned_reads[sample_name + '_forward_strand'], bins =read_length, range=[0,max(aligned_reads['LN'])])
    Freq_fwrd_trim, bins_fwrd_trim = np.histogram(forward_strand[sample_name + '_forward_strand'], bins =read_length, range=[0, 20000])
    #dataframe is created for bins
    hist_fwrd = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd, 'bins_foward': bins_fwrd[:-1]})
    hist_fwrd_trim = pd.DataFrame({'sample_Name': sample_name, 'Freq_forward': Freq_fwrd_trim, 'bins_foward': bins_fwrd_trim[:-1]})
    #Each sample dataframe is concatinated
    fwrd_allsamp_hist = pd.concat([fwrd_allsamp_hist, hist_fwrd])
    fwrd_allsamp_hist_trim = pd.concat([fwrd_allsamp_hist_trim, hist_fwrd_trim])
    fwrd_allsamp_hist.to_csv(snakemake.output['fwrd_allsamp_hist'], sep='\t', index=False)
    fwrd_allsamp_hist_trim.to_csv(snakemake.output['fwrd_allsamp_hist_trim'], sep='\t', index=False)

#Reverse strand
reverse_strand = merge_data.loc[merge_data['Strand'] == '-1']

#Minimum aligned start postion is taken
read_min = reverse_strand.loc[reverse_strand.groupby(['read_Name','reads'])['Start'].idxmin()]

if read_min['Transcript'].str.contains(transcript_ids).any():
    rev_filtered_transcript_data = read_min.query("Transcript == @transcript_ids")
    Freq_rev_fil, bins_rev_fil = np.histogram(read_min['Start'], bins =read_length, range=[0,max(read_min['LN'])])
    hist_rev_fil = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev_fil, 'bins_rev': bins_rev_fil[:-1]})
    rev_allsamp_hist_fil = pd.concat([rev_allsamp_hist_fil, hist_rev_fil])

elif transcript_ids == "all":
    #Values are added to corresponding bins
    Freq_rev, bins_rev = np.histogram(read_min['Start'], bins =read_length, range=[0,max(read_min['LN'])])
    Freq_rev_trim, bins_rev_trim = np.histogram(read_min['Start'], bins =read_length, range=[0, 20000])

    #dataframe is created for bins
    hist_rev = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev, 'bins_rev': bins_rev[:-1]})
    hist_rev_trim = pd.DataFrame({'sample_Name': sample_name, 'Freq_rev': Freq_rev_trim, 'bins_rev': bins_rev_trim[:-1]})

    #Each sample dataframe is concatinated
    rev_allsamp_hist = pd.concat([rev_allsamp_hist, hist_rev])
    rev_allsamp_hist_trim = pd.concat([rev_allsamp_hist_trim, hist_rev_trim])
    rev_allsamp_hist.to_csv(snakemake.output['rev_allsamp_hist'], sep='\t', index=False)
    rev_allsamp_hist_trim.to_csv(snakemake.output['rev_allsamp_hist_trim'], sep='\t', index=False)


if transcript_ids != "all":
    if fwrd_allsamp_hist_fil.empty:
        with open(snakemake.output['fwrd_allsamp_hist_fil'], "w") as fwrd_csvfile:
            print("no aligned transcript",file=fwrd_csvfile)
    else:
        with open(snakemake.output['fwrd_allsamp_hist_fil'], "w") as fwrd_csvfile:
            fwrd_allsamp_hist_fil.to_csv(fwrd_csvfile, sep='\t', index=False)
    if rev_allsamp_hist_fil.empty:
        with open(snakemake.output['rev_allsamp_hist_fil'], "w") as rev_csvfile:
            print("no aligned transcript",file=rev_csvfile)
    else:
        with open(snakemake.output['rev_allsamp_hist_fil'], "w") as rev_csvfile:
            rev_allsamp_hist_fil.to_csv(rev_csvfile, sep='\t', index=False)    

