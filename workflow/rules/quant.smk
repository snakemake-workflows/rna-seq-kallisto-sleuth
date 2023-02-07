rule kallisto_index:
    input:
        fasta="resources/transcriptome.cdna.fasta",
    output:
        index="results/kallisto/transcripts.idx",
    params:
        extra="",  # optional parameters
    log:
        "results/logs/kallisto/index.log",
    wrapper:
        "v1.23.1/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=get_trimmed,
        index="results/kallisto/transcripts.idx",
    output:
        directory("results/kallisto/{sample}-{unit}"),
    log:
        "results/logs/kallisto/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    # around 4 gb of memory usage with 1 thread. (hg38)
    # over 8gb peak memory usage with 8 threads. (hg38)
    threads: 8
    resources:
        mem_mb=10000,
    wrapper:
        "v1.23.1/bio/kallisto/quant"
