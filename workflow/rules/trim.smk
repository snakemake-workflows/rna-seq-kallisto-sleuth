rule cutadapt_pe:
    input:
        get_fastqs
    output:
        fastq1="analysis/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="analysis/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="analysis/trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "analysis/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.31.1/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastqs
    output:
        fastq="analysis/trimmed/{sample}-{unit}.fastq.gz",
        qc="analysis/trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "analysis/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.31.1/bio/cutadapt/se"
