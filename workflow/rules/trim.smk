rule cutadapt_pe:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        "{}".format(config["params"]["cutadapt-pe"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "0.31.1/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastqs,
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        "{}".format(config["params"]["cutadapt-se"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "0.31.1/bio/cutadapt/se"

if config["experiment"]["is-3-prime-rna-seq"]:

    rule cutadapt1:
        input:
            get_fastqs,
        output:
            fastq="results/trimmed/{sample}-{unit}.1.1.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.1.1.qc.txt",
        params:
            extra=lambda w: str( "-m 20 -O 20 -a ""polyA=A{20}"" -a ""QUALITY=G{20}"" -n 2")
        log:
            "results/logs/cutadapt/{sample}-{unit}.1.1.log",
        wrapper:
            "v1.14.1/bio/cutadapt/se"
    
    rule cutadapt2:
        input:
            fastq="results/trimmed/{sample}-{unit}.1.1.fastq.gz",
        output:
            fastq="results/trimmed/{sample}-{unit}.2.1.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.2.1.qc.txt",
        params:
            adapters=lambda w: str("-m 20 -O 3 --nextseq-trim=10 -a \"""A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000\""""),
        log:
            "results/logs/cutadapt/{sample}-{unit}.2.1.log",
        wrapper:
            "v1.14.1/bio/cutadapt/se"
    
    rule cutadapt3:
        input:
            fastq="results/trimmed/{sample}-{unit}.2.1.fastq.gz",
        output:
            fastq="results/trimmed/{sample}-{unit}.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.qc.txt",
        params:
            adapters=lambda w: str("-m 20 -O 20 -g \"""r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20\""" --discard-trimmed"),
        log:
            "results/logs/cutadapt/{sample}-{unit}.log",
        wrapper:
            "v1.14.1/bio/cutadapt/se"


rule get_max_read_length:
    input:
        get_all_fastqs,
    output:
        "results/stats/max-read-length.json",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-max-read-length.py"
        