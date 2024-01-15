import pandas as pd
import numpy as np
import sys
import json


sys.stderr = open(snakemake.log[0], "w")

# Get the read-length
f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()

# Get the transcript IDs if executing for individual transcripts or else for all transcripts
transcript_ids = snakemake.params["each_transcript"]

# Define data-frame for each strand
fwrd_allsamp_hist = pd.DataFrame([])
fwrd_allsamp_hist_trim = pd.DataFrame([])
fwrd_allsamp_hist_fil = pd.DataFrame([])
rev_allsamp_hist = pd.DataFrame([])
rev_allsamp_hist_trim = pd.DataFrame([])
rev_allsamp_hist_fil = pd.DataFrame([])

# Get the sample names
sample_name = snakemake.params["samples"]

# TSV file reading
trans_length_data = pd.read_csv(
    snakemake.input["transcripts_annotation"],
    sep="\t",
).drop(columns=["main_transcript_per_gene"])

# Aligned text file reading
align_bam_txt = pd.read_csv(
    snakemake.input["aligned_file"],
    sep="\t",
    names=["read_Name", "transcript", "start", "reads", "quality"],
)

# Both transcript len and start postion are merged based on same transcript ID
merge_data = align_bam_txt.merge(trans_length_data, on="transcript")

# Forward strand
forward_strand = merge_data.loc[merge_data["strand"] == "+"]

# Each read postion is calcuated
forward_strand[sample_name + "_forward_strand"] = (
    forward_strand["transcript_length"] - forward_strand["start"]
)
aligned_reads = forward_strand.loc[
    forward_strand.groupby(["read_Name", "reads"])[
        sample_name + "_forward_strand"
    ].idxmin()
]
if aligned_reads["transcript"].str.contains(transcript_ids).any():
    # Get aligned read postion of the given transcript
    fwrd_filtered_transcript_data = aligned_reads.query("transcript == @transcript_ids")
    Freq_fwrd_fil, bins_fwrd_fil = np.histogram(
        fwrd_filtered_transcript_data[sample_name + "_forward_strand"],
        bins=read_length,
        range=[0, max(fwrd_filtered_transcript_data["transcript_length"])],
    )
    hist_fwrd_fil = pd.DataFrame(
        {
            "sample_Name": sample_name,
            "Freq_forward": Freq_fwrd_fil,
            "bins_foward": bins_fwrd_fil[:-1],
        }
    )
    fwrd_allsamp_hist_fil = pd.concat([fwrd_allsamp_hist_fil, hist_fwrd_fil])
elif transcript_ids == "all":
    # Values added to corresponding bins
    Freq_fwrd, bins_fwrd = np.histogram(
        aligned_reads[sample_name + "_forward_strand"],
        bins=read_length,
        range=[0, max(aligned_reads["transcript_length"])],
    )
    Freq_fwrd_trim, bins_fwrd_trim = np.histogram(
        forward_strand[sample_name + "_forward_strand"],
        bins=read_length,
        range=[0, 20000],
    )
    # Dataframe created for bins
    hist_fwrd = pd.DataFrame(
        {
            "sample_Name": sample_name,
            "Freq_forward": Freq_fwrd,
            "bins_foward": bins_fwrd[:-1],
        }
    )
    hist_fwrd_trim = pd.DataFrame(
        {
            "sample_Name": sample_name,
            "Freq_forward": Freq_fwrd_trim,
            "bins_foward": bins_fwrd_trim[:-1],
        }
    )
    # Each sample dataframe is concatinated
    fwrd_allsamp_hist = pd.concat([fwrd_allsamp_hist, hist_fwrd])
    fwrd_allsamp_hist_trim = pd.concat([fwrd_allsamp_hist_trim, hist_fwrd_trim])
    fwrd_allsamp_hist.to_csv(
        snakemake.output["fwrd_allsamp_hist"], sep="\t", index=False
    )
    fwrd_allsamp_hist_trim.to_csv(
        snakemake.output["fwrd_allsamp_hist_trim"], sep="\t", index=False
    )

# Reverse strand
reverse_strand = merge_data.loc[merge_data["strand"] == "-"]

# Minimum aligned start postion is taken
read_min = reverse_strand.loc[
    reverse_strand.groupby(["read_Name", "reads"])["start"].idxmin()
]

if read_min["transcript"].str.contains(transcript_ids).any():
    # Get aligned read postion of the given transcript
    rev_filtered_transcript_data = read_min.query("transcript == @transcript_ids")
    Freq_rev_fil, bins_rev_fil = np.histogram(
        rev_filtered_transcript_data["start"],
        bins=read_length,
        range=[0, max(rev_filtered_transcript_data["transcript_length"])],
    )
    hist_rev_fil = pd.DataFrame(
        {
            "sample_Name": sample_name,
            "Freq_rev": Freq_rev_fil,
            "bins_rev": bins_rev_fil[:-1],
        }
    )
    rev_allsamp_hist_fil = pd.concat([rev_allsamp_hist_fil, hist_rev_fil])

elif transcript_ids == "all":
    # Values added to corresponding bins
    Freq_rev, bins_rev = np.histogram(
        read_min["start"],
        bins=read_length,
        range=[0, max(read_min["transcript_length"])],
    )
    Freq_rev_trim, bins_rev_trim = np.histogram(
        read_min["start"], bins=read_length, range=[0, 20000]
    )

    # Dataframe created for bins
    hist_rev = pd.DataFrame(
        {"sample_Name": sample_name, "Freq_rev": Freq_rev, "bins_rev": bins_rev[:-1]}
    )
    hist_rev_trim = pd.DataFrame(
        {
            "sample_Name": sample_name,
            "Freq_rev": Freq_rev_trim,
            "bins_rev": bins_rev_trim[:-1],
        }
    )

    # Each sample dataframe is concatinated
    rev_allsamp_hist = pd.concat([rev_allsamp_hist, hist_rev])
    rev_allsamp_hist_trim = pd.concat([rev_allsamp_hist_trim, hist_rev_trim])
    rev_allsamp_hist.to_csv(snakemake.output["rev_allsamp_hist"], sep="\t", index=False)
    rev_allsamp_hist_trim.to_csv(
        snakemake.output["rev_allsamp_hist_trim"], sep="\t", index=False
    )


# Write bins to each file
if transcript_ids != "all":
    if fwrd_allsamp_hist_fil.empty:
        with open(snakemake.output["fwrd_allsamp_hist_fil"], "w") as fwrd_csvfile:
            print("no reads aligned", file=fwrd_csvfile)
    else:
        with open(snakemake.output["fwrd_allsamp_hist_fil"], "w") as fwrd_csvfile:
            fwrd_allsamp_hist_fil.to_csv(fwrd_csvfile, sep="\t", index=False)
    if rev_allsamp_hist_fil.empty:
        with open(snakemake.output["rev_allsamp_hist_fil"], "w") as rev_csvfile:
            print("no reads aligned", file=rev_csvfile)
    else:
        with open(snakemake.output["rev_allsamp_hist_fil"], "w") as rev_csvfile:
            rev_allsamp_hist_fil.to_csv(rev_csvfile, sep="\t", index=False)
