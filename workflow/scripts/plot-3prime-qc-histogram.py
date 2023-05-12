from turtle import title
import altair as alt
from altair_saver import save
import pandas as pd
import numpy as np
import pysam
from scipy import stats
import sys
import json
import glob
import os

sys.stderr = open(snakemake.log[0], "w")

# reading the file read length
f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()

# Get the transcript IDs if executing for individual transcripts or else for all transcripts
transcript_ids = snakemake.params["each_transcript"]

# Getting the bin files
fwrd_full = []
fwrd_trim = []
fwrd_fil = []

rev_full = []
rev_trim = []
rev_fil = []

# Append all sample bin files into a single dataframe
if transcript_ids != "all":
    fwrd_allsamp_hist_fil = snakemake.input["fwrd_allsamp_hist_fil"]
    rev_allsamp_hist_fil = snakemake.input["rev_allsamp_hist_fil"]
    if fwrd_allsamp_hist_fil != "":
        for filename in fwrd_allsamp_hist_fil:
            df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
            fwrd_fil.append(
                df.loc[
                    :, df.columns.isin(["sample_Name", "Freq_forward", "bins_foward"])
                ]
            )
        fwrd_allsamp_hist_fil = pd.concat(fwrd_fil, axis=0, ignore_index=True)

    if rev_allsamp_hist_fil != "":
        for filename in rev_allsamp_hist_fil:
            df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
            rev_fil.append(
                df.loc[:, df.columns.isin(["sample_Name", "Freq_rev", "bins_rev"])]
            )
        rev_allsamp_hist_fil = pd.concat(rev_fil, axis=0, ignore_index=True)
else:
    file_fwrd_allsamp_hist = snakemake.input["fwrd_allsamp_hist"]
    file_fwrd_allsamp_hist_trim = snakemake.input["fwrd_allsamp_hist_trim"]
    file_rev_allsamp_hist = snakemake.input["rev_allsamp_hist"]
    file_rev_allsamp_hist_trim = snakemake.input["rev_allsamp_hist_trim"]

    for filename in file_fwrd_allsamp_hist:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        fwrd_full.append(df)
    fwrd_allsamp_hist = pd.concat(fwrd_full, axis=0, ignore_index=True)

    for filename in file_fwrd_allsamp_hist_trim:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        fwrd_trim.append(df)
    fwrd_allsamp_hist_trim = pd.concat(fwrd_trim, axis=0, ignore_index=True)

    for filename in file_rev_allsamp_hist:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        rev_full.append(df)
    rev_allsamp_hist = pd.concat(rev_full, axis=0, ignore_index=True)

    for filename in file_rev_allsamp_hist_trim:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        rev_trim.append(df)
    rev_allsamp_hist_trim = pd.concat(rev_trim, axis=0, ignore_index=True)

# Plot the qc-histogram

if transcript_ids != "all":
    if not fwrd_allsamp_hist_fil.empty:

        hist_fwrd_full = (
            alt.Chart(fwrd_allsamp_hist_fil)
            .mark_line(interpolate="step-after")
            .encode(
                x=alt.X(
                    "bins_foward",
                    title="difference between transcript length and read start",
                ),
                y=alt.Y("Freq_forward:Q", title="Count of Records"),
                color="sample_Name",
            )
            .properties(title="forward strand transcripts (full length)")
        )

        fwrd_allsamp_hist_fil["read_length"] = read_length

        cht_rd_len_fwd_full = (
            alt.Chart(fwrd_allsamp_hist_fil)
            .mark_rule(color="black", strokeDash=[3, 5])
            .encode(x="read_length")
        )

        Fwd_chart_full = hist_fwrd_full + cht_rd_len_fwd_full
        Fwd_chart_full.save(snakemake.output["full_sample_QC"])

    elif not rev_allsamp_hist_fil.empty:
        hist_rev_full = (
            alt.Chart(rev_allsamp_hist_fil)
            .mark_line(interpolate="step-after")
            .encode(
                x=alt.X("bins_rev", title="read start position in the transcript"),
                y=alt.Y("Freq_rev:Q", title="Count of Records"),
                color="sample_Name",
            )
            .properties(title="Reverse strand transcripts (full length)")
        )

        rev_allsamp_hist_fil["read_length"] = read_length

        cht_rd_len_rev = (
            alt.Chart(rev_allsamp_hist_fil)
            .mark_rule(color="black", strokeDash=[3, 5])
            .encode(x="read_length")
        )

        Rev_chart_full = hist_rev_full + cht_rd_len_rev

        Rev_chart_full.save(snakemake.output["full_sample_QC"])

    elif fwrd_allsamp_hist_fil.empty:
        # Configure empty histogram (forward strand)
        empty_histogram = alt.Chart().mark_text(text="no reads aligned", size=20)
        empty_histogram.save(snakemake.output["full_sample_QC"])
    # Configure empty histogram (reverse strand)
    elif rev_allsamp_hist_fil.empty:
        empty_histogram = alt.Chart().mark_text(text="no reads aligned", size=20)
        empty_histogram.save(snakemake.output["full_sample_QC"])


else:
    # Histogram for forward strand
    # Histogram for full len transcript
    hist_fwrd_full = (
        alt.Chart(fwrd_allsamp_hist)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X(
                "bins_foward",
                title="difference between transcript length and read start",
            ),
            y=alt.Y("Freq_forward:Q", title="Count of Records"),
            color="sample_Name",
        )
        .properties(title="forward strand transcripts (full length)")
    )
    fwrd_allsamp_hist["read_length"] = read_length
    cht_rd_len_fwd_full = (
        alt.Chart(fwrd_allsamp_hist)
        .mark_rule(color="black", strokeDash=[3, 5])
        .encode(x="read_length")
    )
    Fwd_chart_full = hist_fwrd_full + cht_rd_len_fwd_full

    # Histogram plot for 20000 bp len transcript
    hist_fwrd_trimd = (
        alt.Chart(fwrd_allsamp_hist_trim)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X(
                "bins_foward",
                title="difference between transcript length and read start",
            ),
            y=alt.Y("Freq_forward:Q", title="Count of Records"),
            color="sample_Name",
        )
        .properties(title="forward strand transcripts (showing 1-20000bp)")
    )
    fwrd_allsamp_hist_trim["read_length"] = read_length
    cht_rd_len_fwd_trim = (
        alt.Chart(fwrd_allsamp_hist_trim)
        .mark_rule(color="black", strokeDash=[3, 5])
        .encode(x="read_length")
    )
    Fwd_chart_trim = hist_fwrd_trimd + cht_rd_len_fwd_trim

    Fwd_chart = alt.hconcat(Fwd_chart_trim, Fwd_chart_full)

    # Histogram for reverse strand
    # Histogram for full len transcript
    hist_rev_full = (
        alt.Chart(rev_allsamp_hist)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("bins_rev", title="read start position in the transcript"),
            y=alt.Y("Freq_rev:Q", title="Count of Records"),
            color="sample_Name",
        )
        .properties(title="Reverse strand transcripts (full length)")
    )
    rev_allsamp_hist["read_length"] = read_length
    cht_rd_len_rev = (
        alt.Chart(rev_allsamp_hist)
        .mark_rule(color="black", strokeDash=[3, 5])
        .encode(x="read_length")
    )
    rev_chart_full = hist_rev_full + cht_rd_len_rev

    # Histogram plot for 20000 bp len transcript
    hist_rev_trim = (
        alt.Chart(rev_allsamp_hist_trim)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("bins_rev", title="read start position in the transcript"),
            y=alt.Y("Freq_rev:Q", title="Count of Records"),
            color="sample_Name",
        )
        .properties(title="Reverse strand transcripts (showing 1-20000bp)")
    )
    rev_allsamp_hist_trim["read_length"] = read_length
    cht_rd_len_rev_trim = (
        alt.Chart(rev_allsamp_hist_trim)
        .mark_rule(color="black", strokeDash=[3, 5])
        .encode(x="read_length")
    )
    rev_chart_trim = hist_rev_trim + cht_rd_len_rev_trim

    Rev_chart = alt.hconcat(rev_chart_trim, rev_chart_full)

    Final_chart = alt.vconcat(Fwd_chart, Rev_chart)
    Final_chart.save(snakemake.output["full_sample_QC"])
