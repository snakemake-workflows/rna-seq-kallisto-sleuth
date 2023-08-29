import sys
import re
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
        seq_record.seq = re.sub("TTTT+$|AAAA+$", "", str(seq_record.seq))
        SeqIO.write(
            seq_record,
            transcript_clean_cdna_fasta,
            "fasta"
        )
