rule cutadapt_pe:
    input:
        get_fastqs
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.31.1/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastqs
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.31.1/bio/cutadapt/se"
