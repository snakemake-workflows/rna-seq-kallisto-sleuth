if is_3prime_experiment and three_prime_vendor == "lexogen":

    # Rule cutadapt1 checks and removes poly-A tails and sequence qualilty score <20.
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `-O 20 -a ""polyA=A{20}"" -a ""QUALITY=G{20}"" -n 2`:
    #   * Removes reads having polyA (18bp) or polyG (18bp) and repeat the process twice by -n 2
    # In short, it removes reads that are predominantly polyA or polyG.

    # Rule cutadapt2 checks if reads contains the polyA-stretch + adapter at the 3'
    # end with minimum overlap 3bp, maximum error of 1bp every 10bp
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `--nextseq-trim=10`
    #   * Quality score<20 is trimmed, based on NextSeq/NovaSeq 2-color SBS system
    #   * `-a A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000`

    # Rule cutadapt3 removes final set of adapters from 5'end and discards reads based on trimmed read length
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `-O 20 -g r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20 --discard-trimmed"`
    #   * Trim the 5' end of adapter with miminum overlap of 20bp and discard its read by --discard-trimmed

    # https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-1
    rule cutadapt1:
        input:
            get_fastqs,
        output:
            fastq="results/trimmed/{sample}-{unit}.1.1.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.1.1.qc.txt",
        params:
            extra=lambda w: str(
                "-m 20 -O 20 -a " "polyA=A{20}" " -a " "QUALITY=G{20}" " -n 2"
            ),
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
            adapters=lambda w: str(
                '-m 20 -O 3 --nextseq-trim=10 -a "'
                'A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"'
                ""
            ),
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
            adapters=lambda w: str(
                '-m 20 -O 20 -g "'
                'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"'
                " --discard-trimmed"
            ),
        log:
            "results/logs/cutadapt/{sample}-{unit}.log",
        wrapper:
            "v1.14.1/bio/cutadapt/se"


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
