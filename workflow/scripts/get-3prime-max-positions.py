import pandas as pd
import numpy as np
import pysam
import sys

sys.stderr = open(snakemake.log[0], "w")

bam_file = snakemake.input["canonical_mapped_bam"]
sample_name = bam_file.split(".canonical.mapped.sorted.bam")[0]
# Bam file reading
bam_file = pysam.AlignmentFile(snakemake.input["canonical_mapped_bam"], "rb")
bam_header = bam_file.header.to_dict()
trans_length_data = pd.DataFrame(bam_header.get("SQ"))
trans_length_data.rename(columns={"SN": "Transcript_ID"}, inplace=True)

# Aligned text file reading
align_bam_txt = pd.read_csv(
    snakemake.input["canonical_mapped_pos"],
    sep="\t",
    names=["read_name", "Transcript_ID", "Start", "read", "Quality"],
)
align_bam_txt["Strand"] = align_bam_txt["Transcript_ID"].str.split("_", 1).str[1]
align_bam_txt["Transcript"] = align_bam_txt["Transcript_ID"].str.split("_", 1).str[0]
merge_data = align_bam_txt.merge(trans_length_data, on="Transcript_ID")

# reads aligned to forward strand
forward_strand = merge_data.loc[merge_data["Strand"] == "1"]
forward_strand[sample_name + "_forward_strand"] = (
    forward_strand["LN"] - forward_strand["Start"]
)
aligned_reads = forward_strand.loc[
    forward_strand.groupby(["read_name", "read"])[
        sample_name + "_forward_strand"
    ].idxmin()
]

# reads aligned to reverse strand
reverse_strand = merge_data.loc[merge_data["Strand"] == "-1"]
read_min = reverse_strand.loc[
    reverse_strand.groupby(["read_name", "read"])["Start"].idxmin()
]
aligned_reads = pd.concat([aligned_reads, read_min])
aligned_reads.to_csv(
    snakemake.output["canonical_mapped_3prime_pos"],
    columns=["read_name"],
    sep="\t",
    index=False,
)
