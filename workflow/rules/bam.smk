rule bam_paired_to_fastq:
    input:
        lookup(query="sample == '{sample}' & unit == '{unit}'", within=units, cols="bam_paired")
    output:
        "results/{sample}-{unit}.1.fq.gz",
        "results/{sample}-{unit}.1.fq.gz",
    log:
        "logs/{sample}-{unit}.separate.log",
    params:
        fastq="-n",
    threads: 3
    wrapper:
        "v3.10.2/bio/samtools/fastq/separate"


rule bam_single_to_fastq:
    input:
        lookup(query="sample == '{sample}' & unit == '{unit}'", within=units, cols="bam_single")
    output:
        "results/{sample}-{unit}.fq.gz",
    log:
        "logs/{sample}-{unit}.interleaved.log",
    threads: 3
    wrapper:
        "v3.10.2/bio/samtools/fastq/interleaved"
