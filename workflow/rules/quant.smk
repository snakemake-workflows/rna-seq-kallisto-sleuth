rule kallisto_cds_index:
    input:
        fasta="resources/transcriptome_clean.cdna.fasta"
            if config["experiment"]["3-prime-rna-seq"]["activate"]
            else "resources/transcriptome.cdna.fasta",                
    output:
        index="results/kallisto_cds/transcripts.idx",
    log:
        "results/logs/kallisto_cds/index.log",
    threads: 1
    wrapper:
        "v1.17.4/bio/kallisto/index"


rule kallisto_cds_quant:
    input:
        fastq=get_trimmed,
        index="results/kallisto_cds/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_cds/{sample}-{unit}"),
    log:
        "results/logs/kallisto_cds/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    wrapper:
        "v1.17.4/bio/kallisto/quant"


rule bwa_index:
    input:
        "resources/transcriptome_clean.cdna.fasta"
        if config["experiment"]["3-prime-rna-seq"]["activate"]
        else "resources/transcriptome.cdna.fasta", 
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
        "results/mapped_mem/{sample}-{unit}.bam",
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.17.2/bio/bwa/mem"


rule kallisto_3prime_index:
    input:
        fasta="resources/transcriptome_canonical.3prime.fasta",
    output:
        index="results/kallisto_3prime/transcripts.idx",
    log:
        "results/logs/kallisto_3prime/index.log",
    threads: 1
    wrapper:
        "v1.17.4/bio/kallisto/index"


rule kallisto_3prime_quant:
    input:
        fastq="results/canonical_reads/{sample}-{unit}.fastq",
        index="results/kallisto_3prime/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_3prime/{sample}-{unit}"),
    log:
        "results/logs/kallisto_3prime/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    wrapper:
        "v1.17.4/bio/kallisto/quant"


rule mapped_bam:
    input:
        bam_file="results/mapped_mem/{sample}-{unit}.bam",
    output:
        mapped_bam=temp("results/mapped_bam/{sample}-{unit}.mapped.bam"),
    log:
        "results/logs/mapped_bam/{sample}-{unit}.mapped-bam.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b -F 4 {input.bam_file} > {output.mapped_bam} 2> {log}"


rule get_mapped_canonical_transcripts:
    input:
        mapped_bam="results/mapped_bam/{sample}-{unit}.mapped.bam",
        canonical_ids="resources/canonical_ids.csv",
    output:
        canonical_mapped_bam=temp("results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam"),
    log:
        "results/logs/canonical_mapped_bam/{sample}-{unit}.canonical-mapped-bam.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -h {input.mapped_bam} |  cut -f1-12 | grep -f {input.canonical_ids} > {output.canonical_mapped_bam}  2> {log}"


rule get_mapped_canonical_positions:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
    output:
        canonical_mapped_pos=temp("results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.position.txt"),
    log:
        "results/logs/canonical_mapped_bam/{sample}-{unit}.canonical-mapped-pos.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {input.canonical_mapped_bam} | cut -f1,3,4,10,11  > {output}  2> {log}"


rule get_closest_3prime_aligned_pos:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
        canonical_mapped_pos="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.position.txt",
    output:
        canonical_mapped_3prime_pos=temp("results/mapped_3prime_bam/{sample}-{unit}.canonical.mapped.3prime_pos.txt"),
        bam_sorted=temp("results/mapped_3prime_bam/{sample}-{unit}.canonical.sorted.mapped.pos.bam")
    params:
        samples="results/kallisto_cds/{sample}-{unit}",
    log:
        "results/logs/mapped_3prime_bam/{sample}-{unit}.mapped.pos.log",
    conda:
        "../envs/QC.yaml",
    script:
        "../scripts/get-3prime-max-positions.py"


rule get_closest_3prime_aligned_pos_bam:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
        canonical_mapped_3prime_pos="results/mapped_3prime_bam/{sample}-{unit}.canonical.mapped.3prime_pos.txt",
    output:
        canonical_mapped_3prime_bam=temp("results/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.bam"),
    log:
        "results/logs/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -R {input.canonical_mapped_3prime_pos} {input.canonical_mapped_bam} -o {output.canonical_mapped_3prime_bam}  2> {log}"


rule get_canonical_fastq:
    input:
        canonical_3prime_mapped_bam="results/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.bam",
    output:
        canonical_fastq="results/canonical_reads/{sample}-{unit}.fastq",
    log:
        "results/logs/canonical_3prime_mapped_bam/{sample}-{unit}.canonical_3prime_mapped.fastq.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools bam2fq {input.canonical_3prime_mapped_bam} > {output.canonical_fastq}  2> {log}"


rule get_aligned_pos:
    input:
        bam_file="results/kallisto_cds/{sample}-{unit}",
    output:
        aligned_files=temp("results/QC/{sample}-{unit}.aligned.txt"),
    log:
        "results/logs/QC/{sample}-{unit}.aligned.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {input.bam_file}/pseudoalignments.bam | cut -f1,3,4,10,11  > {output} 2> {log}"  


rule get_read_dist:
    input:
        aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples()),
        read_length="results/stats/max-read-length.json",
    output:
        report(
            "results/plots/QC/{model}.QC-plot.html",
            caption="../report/plot-QC.rst",
            category="QC",
        ),
    params:
        samples=expand("results/kallisto_cds/{unit.sample}-{unit.unit}",unit=units.itertuples()),
    log:
        "results/logs/QC/{model}.QC-plot.log",
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/plot-sample-QC-histogram.py"


if config["experiment"]["3-prime-rna-seq"]["plot-qc"] != "all":
    
    
    rule get_ind_transcript_histograms:
        input:
            aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples()),
        output:
            report(
                "results/QC/{ind_transcripts}.QC_plot.html",
                caption="../report/plot-QC.rst",
                category="QC",
            ),
        params:
            each_transcript = "{ind_transcripts}",
            read_length="results/stats/max-read-length.json",
            #samples=expand("results/kallisto_cds/{unit.sample}-{unit.unit}",unit=units.itertuples()),
            samples=expand("results/kallisto_cds/Mel-86c_p49_Ctrl-1",unit=units.itertuples()),
        log:
            "results/logs/QC/{ind_transcripts}.QC_plot.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/plot_ind-transcripts_histogram.py"