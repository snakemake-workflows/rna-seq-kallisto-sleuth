import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    all_transcripts = pd.read_csv(
        snakemake.input["main_transcripts_per_gene"],
        sep="\t",
    )
    main_transcripts_per_gene = set(
        all_transcripts.loc[
            all_transcripts["main_transcript_per_gene"] == 1, "transcript"
        ]
    )

    for seq_record in SeqIO.parse(snakemake.input["fasta"], "fasta"):
        if seq_record.id in main_transcripts_per_gene:
            SeqIO.write(seq_record, transcript_clean_cdna_fasta, "fasta")
