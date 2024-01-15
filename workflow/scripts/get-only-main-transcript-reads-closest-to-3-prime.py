import sys

import pandas as pd
import pysam

sys.stderr = open(snakemake.log[0], "w")

transcript_annotations = pd.read_csv(
    snakemake.input["annotation"],
    sep="\t",
)


def print_read_closest_to_3_prime(read_set: set[pysam.AlignedSegment]) -> None:
    min_dist = -1
    for read in read_set:
        if read.is_mapped:
            transcript = transcript_annotations.loc[
                transcript_annotations["transcript"] == read.reference_name,
                ["transcript_length", "main_transcript_per_gene"],
            ].reset_index()
            if len(transcript) == 0:
                sys.stderr.write(
                    f"Warning: Transcript '{read.reference_name}' not found in downloaded annotations. Skipping read '{read.query_name}' that maps to it."
                )
                continue
            # use read start distance, as reference skips increase that distance
            # and also indicate a suboptimal alignment
            distance = transcript.at[0, "transcript_length"] - read.reference_start
            if min_dist == -1 or distance < min_dist:
                min_dist = distance
                min_dist_read = read
                min_dist_transcript = transcript
    if min_dist != -1 and min_dist_transcript.at[0, "main_transcript_per_gene"] == 1:
        records_out.write(min_dist_read)


with pysam.AlignmentFile(
    snakemake.input["bam"], "rb"
) as records_in, pysam.AlignmentFile(
    snakemake.output["main_transcript_reads_closest_to_3_prime"],
    "wb",
    template=records_in,
) as records_out:
    record_iterator = records_in.fetch(until_eof=True)
    current_read = next(record_iterator)
    current_queryname_set = set([current_read])
    while (next_read := next(record_iterator, None)) is not None:
        if next_read.query_name != current_read.query_name:
            print_read_closest_to_3_prime(current_queryname_set)
            current_queryname_set = set([next_read])
        else:
            current_queryname_set.add(next_read)
        current_read = next_read
    print_read_closest_to_3_prime(current_queryname_set)
