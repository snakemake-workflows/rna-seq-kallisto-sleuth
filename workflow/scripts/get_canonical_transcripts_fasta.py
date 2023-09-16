import sys
import re
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    mane_select_transcripts = set(
        pd.read_csv(
            snakemake.input["mane_select_transcripts"],
            sep="\t",
            usecols=["transcript"],
        )
    )

    for seq_record in SeqIO.parse(snakemake.input["fasta"], "fasta"):
        if seq_record.id in mane_select_transcripts:
            SeqIO.write(
                seq_record,
                transcript_clean_cdna_fasta,
                "fasta"
            )
