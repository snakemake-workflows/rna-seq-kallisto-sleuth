import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    all_transcripts = pd.read_csv(
            snakemake.input["mane_select_transcripts"],
            sep="\t",
        )
    mane_select_transcripts = set(
        all_transcripts.loc[
            all_transcripts["transcript_mane_select"] == 1,
            "transcript"
        ] 
    )

    for seq_record in SeqIO.parse(snakemake.input["fasta"], "fasta"):
        if seq_record.id in mane_select_transcripts:
            SeqIO.write(
                seq_record,
                transcript_clean_cdna_fasta,
                "fasta"
            )
