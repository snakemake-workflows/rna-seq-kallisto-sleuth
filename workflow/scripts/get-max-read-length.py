import json
import pysam

max_read_length = max(
    len(rec.sequence) for f in snakemake.input for rec in pysam.FastxFile(f)
)

with open(snakemake.output[0], "w") as out:
    json.dump(max_read_length, out)
