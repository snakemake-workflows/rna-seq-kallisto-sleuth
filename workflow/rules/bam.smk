rule bam_paired_to_fastq:
    input:
        lookup(
            query="sample == '{sample}' & unit == '{unit}'",
            within=units,
            cols="bam_paired",
        ),
    output:
        "results/fastq/{sample}-{unit}.1.fq.gz",
        "results/fastq/{sample}-{unit}.2.fq.gz",
    log:
        "logs/fastq/{sample}-{unit}.separate.log",
    params:
        fastq="-n",
    threads: 3
    wrapper:
        "v3.12.2/bio/samtools/fastq/separate"


rule bam_single_to_fastq:
    input:
        lookup(
            query="sample == '{sample}' & unit == '{unit}'",
            within=units,
            cols="bam_single",
        ),
    output:
        "results/fastq/{sample}-{unit}.fq.gz",
    log:
        "logs/fastq/{sample}-{unit}.interleaved.log",
    threads: 3
    wrapper:
        "v3.12.2/bio/samtools/fastq/interleaved"
