rule kallisto_cds_index:
    input:
<<<<<<< HEAD
        "resources/transcriptome_clean.cdna.fasta"
        if config["experiment"]["is-3-prime-rna-seq"]
        else "resources/transcriptome.cdna.fasta",                
=======
        "resources/transcriptome_clean.cdna.fasta",
>>>>>>> a1cce49705730f2c73aa5869fee0fa86fab6ff7c
    output:
        "results/kallisto_cds/transcripts.idx",
    log:
        "results/logs/kallisto_cds/index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"
    
rule kallisto_cds_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto_cds/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_cds/{sample}-{unit}"),
        bam_file="results/kallisto_cds/{sample}-{unit}/pseudoalignments.bam",
    log:
        "results/logs/kallisto_cds/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"

if config["experiment"]["is-3-prime-rna-seq"]:
    rule kallisto_3prime_index:
        input:
            "resources/transcriptome.3prime.fasta",
        output:
            "results/kallisto_3prime/transcripts.idx",
        log:
            "results/logs/kallisto_3prime/index.log",
        conda:
            "../envs/kallisto.yaml"
        shell:
            "kallisto index -i {output} {input} 2> {log}"

<<<<<<< HEAD



    rule kallisto_3prime_quant:
        input:
            fq=get_trimmed,
            idx="results/kallisto_3prime/transcripts.idx",
        output:
            kallisto_folder=directory("results/kallisto_3prime/{sample}-{unit}"),
        log:
            "results/logs/kallisto_3prime/quant/{sample}-{unit}.log",
        params:
            extra=kallisto_params,
        conda:
            "../envs/kallisto.yaml"
        shell:
            "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"

if config["experiment"]["is-3-prime-rna-seq"]:
    rule get_aligned_pos:
        input:
            bam_file="results/kallisto_cds/{sample}-{unit}/pseudoalignments.bam",
        output:
            aligned_files="results/QC/{sample}-{unit}.aligned.txt",
        log:
            "results/logs/QC/{sample}-{unit}.aligned.log",
        conda:
            "../envs/aligned_pos.yaml"
        shell:
            "samtools view {input.bam_file} | cut -f1,3,4  > {output} 2> {log}"
=======
rule kallisto_3prime_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto_3prime/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_3prime/{sample}-{unit}"),
    log:
        "results/logs/kallisto_3prime/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"


rule get_aligned_pos:
    input:
        bam_file="results/kallisto_cds/{sample}-{unit}/pseudoalignments.bam",
        #expand("results/kallisto/{unit.sample}-{unit.unit}/pseudoalignments.bam", unit=units.itertuples()),
    output:
        aligned_files="results/QC/{sample}-{unit}.aligned.txt",
        #"results/QC/{sample}-{unit}.aligned.txt"
    log:
        "results/logs/QC/{sample}-{unit}.aligned.log",
    shell:
        "samtools view {input.bam_file} | cut -f1,3,4  > {output} 2> {log}"
>>>>>>> a1cce49705730f2c73aa5869fee0fa86fab6ff7c


rule get_read_dist:
    input:
<<<<<<< HEAD
        aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples()),
        read_length="results/stats/max-read-length.json",
    output:
        report(
            "results/QC/QC_plot.html",
            caption="../report/plot-QC.rst",
            category="QC",
        ),
    params:
        samples=expand("results/kallisto_cds/{unit.sample}-{unit.unit}",unit=units.itertuples()),  
=======
        aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples())
        #aligned_file="results/QC/{sample}-{unit}.aligned.txt",
        #bam_path="results/kallisto_cds/{sample}-{unit}/"
    output:
        histogram="results/QC/QC_plot.html",
    params:
        samples=expand("results/kallisto_cds/{unit.sample}-{unit.unit}",unit=units.itertuples()),
        read_length="results/stats/max-read-length.json",
>>>>>>> a1cce49705730f2c73aa5869fee0fa86fab6ff7c
    log:
        "results/logs/QC/QC_plot.log",
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/get_histogram.py"