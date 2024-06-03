rule bam_paired_to_fastq:
    input:
        "{bam_file}.bam",
    output:
        "{bam_file}.1.fq",
        "{bam_file}.2.fq",
    log:
        "{bam_file}.separate.log",
    params:
        fastq="-n",
    threads: 3
    wrapper:
        "v3.10.2/bio/samtools/fastq/separate"


rule samtools_fastq_interleaved:
    input:
        "{bam_file}.bam",
    output:
        "{bam_file}.fq",
    log:
        "{bam_file}.interleaved.log",
    params:
        " ",
    threads: 3
    wrapper:
        "v3.10.2/bio/samtools/fastq/interleaved"
