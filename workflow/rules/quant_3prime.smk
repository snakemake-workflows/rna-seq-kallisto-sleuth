rule bwa_index:
    input:
        "resources/transcriptome.cdna.without_poly_a.fasta",
    output:
        idx=multiext("resources/transcriptome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/transcriptome.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.17.2/bio/bwa/index"


rule bwa_mem:
    input:
        reads=get_trimmed,
        idx=multiext("resources/transcriptome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/mapped_mem/{sample}-{unit}.namesorted.bam",
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.17.2/bio/bwa/mem"


rule get_only_mane_select_reads_closest_to_3_prime:
    input:
        bam="results/mapped_mem/{sample}-{unit}.namesorted.bam",
        annotation="resources/transcripts_annotation.mane_strand_length.tsv",
    output:
        mane_select_reads_closest_to_3_prime=temp(
            "results/mapped_3prime_mane/{sample}-{unit}.mane_select_closest_to_3_prime.bam"
        ),
    log:
        "logs/mapped_3prime_bam/{sample}-{unit}.mapped.pos.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-only-mane-select-reads-closest-to-3-prime.py"


rule get_mane_fastq:
    input:
        bam="results/mapped_3prime_mane/{sample}-{unit}.mane_select_closest_to_3_prime.bam"
    output:
        fastq="results/mane_3prime_reads/{sample}-{unit}.fastq",
    log:
        "logs/mane_3prime_reads/{sample}-{unit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools bam2fq {input.bam} > {output.fastq}  2> {log}"


rule kallisto_3prime_index:
    input:
        fasta="resources/transcriptome.cdna.without_poly_a.mane.fasta",
    output:
        index="results/kallisto_3prime/transcripts.3prime.idx",
    log:
        "logs/kallisto_3prime/index.3prime.log",
    threads: 1
    wrapper:
        "v1.23.1/bio/kallisto/index"


rule kallisto_3prime_quant:
    input:
        fastq=kallisto_quant_input,
        index="results/kallisto_3prime/transcripts.3prime.idx",
    output:
        kallisto_folder=directory("results/kallisto_3prime/{sample}-{unit}"),
    log:
        "logs/kallisto_3prime/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    wrapper:
        "v1.23.1/bio/kallisto/quant"


rule kallisto_samtools_sort:
    input:
        "results/kallisto_cdna/{sample}-{unit}",
    output:
        temp("results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam"),
    log:
        "logs/QC/{sample}-{unit}.sorted.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input}/pseudoalignments.bam > {output} 2> {log}"


rule kallisto_samtools_index:
    input:
        "results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam",
    output:
        temp(
            "results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam.bai"
        ),
    log:
        "logs/QC/{sample}-{unit}.sorted.index.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"