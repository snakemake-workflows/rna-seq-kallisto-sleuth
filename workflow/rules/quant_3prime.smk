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


rule get_only_main_transcript_reads_closest_to_3_prime:
    input:
        bam="results/mapped_mem/{sample}-{unit}.namesorted.bam",
        annotation="resources/transcripts_annotation.main_transcript_strand_length.tsv",
    output:
        sam=temp(
            "results/mapped_3prime_main_transcript/{sample}-{unit}.main_transcript_closest_to_3_prime.sam"
        ),
    log:
        "logs/mapped_3prime_bam/{sample}-{unit}.mapped.pos.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "( "
        " samtools view -H {input.bam} >{output.sam}; "
        " samtools view --exclude-flags unmap {input.bam} | "
        "   awk ' "
        "     BEGIN {{ dist = -1; main = 0; }} "
        "     NR==1 {{ "
        "       for (i=1; i<=NF; i++) {{ "
        "         f[$i] = i "
        "       }} "
        "     }} "
        "     FNR==NR {{ "
        '       t[$(f["transcript"]),"main"] = $(f["main_transcript_per_gene"]); '
        '       t[$(f["transcript"]),"len"] = $(f["transcript_length"]); '
        "       next "
        "     }} "
        "     {{ "
        "       if ($1 != read_id) {{ "
        "         if (main == 1) {{ print read; }}; "
        "         dist = -1; main = 0; "
        "       }}; "
        '       read_id = $1; new_dist = t[$3,"len"] - $4; '
        "       if ( new_dist >= -1 && (dist == -1 || dist > new_dist)) {{ "
        '         dist = new_dist; main = t[$3,"main"]; read = $0; '
        "       }} "
        "     }} "
        "     END {{ "
        "       if (main == 1) {{ print read; }} "
        "     }} "
        "     ' "
        "     {input.annotation} - >> {output.sam} "
        ") 2>{log}"


rule get_main_transcript_fastq:
    input:
        sam="results/mapped_3prime_main_transcript/{sample}-{unit}.main_transcript_closest_to_3_prime.sam",
    output:
        fastq="results/main_transcript_3prime_reads/{sample}-{unit}.fastq",
    log:
        "logs/main_transcript_3prime_reads/{sample}-{unit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq {input.sam} > {output.fastq}  2> {log}"


rule kallisto_3prime_index:
    input:
        fasta="resources/transcriptome.cdna.without_poly_a.main_transcript.fasta",
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
