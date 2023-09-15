import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

sample_name = f"{snakemake.wildcards['sample']}-{snakemake.wildcards['unit']}"

# BED file reading
trans_length_data = pd.read_csv(
    snakemake.input["canonical_ids"],
    sep="\t",
    names=["transcript", "transcript_start", "transcript_length", "name", "transcript_mane_select" "strand"],
).drop(columns = ["transcript_start", "name"])

# Aligned text file reading
align_bam_txt = pd.read_csv(
    snakemake.input["canonical_mapped_pos"],
    sep="\t",
    names=["read_name", "transcript", "start", "read", "quality"],
)
merge_data = align_bam_txt.merge(trans_length_data, on="transcript")

# reads aligned to forward strand
forward_strand = merge_data.loc[merge_data["strand"] == "+"]
forward_strand[sample_name + "_forward_strand"] = (
    forward_strand["transcript_length"] - forward_strand["start"]
)
forward_aligned = forward_strand.loc[
    forward_strand.groupby(["read_name", "read"])[
        sample_name + "_forward_strand"
    ].idxmin()
]

# reads aligned to reverse strand
reverse_strand = merge_data.loc[merge_data["strand"] == "-"]
reverse_aligned = reverse_strand.loc[
    reverse_strand.groupby(["read_name", "read"])["start"].idxmin()
]
aligned_reads = pd.concat([forward_aligned, reverse_aligned])
aligned_reads.to_csv(
    snakemake.output["canonical_mapped_3prime_pos"],
    columns=["read_name"],
    sep="\t",
    index=False,
)
