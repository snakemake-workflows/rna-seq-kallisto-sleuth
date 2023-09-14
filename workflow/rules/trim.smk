rule cutadapt_pe:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    threads: 8
    params:
        adapters=config["params"]["cutadapt-pe"]["adapters"],
        extra=config["params"]["cutadapt-pe"]["extra"],
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "v2.6.0/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastqs,
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    threads: 8
    params:
        adapters=config["params"]["cutadapt-se"]["adapters"],
        extra=config["params"]["cutadapt-se"]["extra"],
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "v2.6.0/bio/cutadapt/se"


rule max_read_length:
    input:
        get_all_fastqs,
    output:
        "results/stats/max-read-length.json",
    log:
        "logs/max-read-length.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-max-read-length.py"
