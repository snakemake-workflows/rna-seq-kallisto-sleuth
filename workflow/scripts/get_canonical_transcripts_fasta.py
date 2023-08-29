import sys
import re
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    canonical_ids = set(
        pd.read_csv(
            snakemake.input["canonical_ids"],
            sep="\t",
            names=["transcript", "transcript_start", "transcript_length", "strand"],
        ).drop(columns = ["transcript_start", "transcript_length", "strand"])
        .loc[:, "transcript"]
    )

    for seq_record in SeqIO.parse(snakemake.input["fasta"], "fasta"):
        if seq_record.id in canonical_ids:
            SeqIO.write(
                seq_record,
                transcript_clean_cdna_fasta,
                "fasta"
            )
