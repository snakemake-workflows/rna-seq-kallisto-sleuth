rule cutadapt_pe:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        extra="{}".format(config["params"]["cutadapt-pe"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "v1.22.0/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastqs,
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        extra="{}".format(config["params"]["cutadapt-se"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "v1.22.0/bio/cutadapt/se"
