import json
import pysam
import sys

sys.stderr = open(snakemake.log[0], "w")

max_read_length = max(
    len(rec.sequence) for f in snakemake.input for rec in pysam.FastxFile(f)
)

with open(snakemake.output[0], "w") as out:
    json.dump(max_read_length, out)
