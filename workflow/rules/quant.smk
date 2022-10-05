rule kallisto_cds_index:
    input:
        "resources/transcriptome_clean.cdna.fasta"
        if config["experiment"]["3-prime-rna-seq"]["activate"]
        else "resources/transcriptome.cdna.fasta",                
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

if config["experiment"]["3-prime-rna-seq"]["activate"]:
    
    rule kallisto_3prime_index:
        input:
            "resources/transcriptome_canonical.3prime.fasta",
        output:
            "results/kallisto_3prime/transcripts.idx",
        log:
            "results/logs/kallisto_3prime/index.log",
        conda:
            "../envs/kallisto.yaml"
        shell:
            "kallisto index -i {output} {input} 2> {log}"

    rule kallisto_3prime_quant:
        input:
            fq="results/canonical_reads/{sample}-{unit}.fastq",
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

if config["experiment"]["3-prime-rna-seq"]["activate"]:
    
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
            "samtools view {input.bam_file} | cut -f1,3,4,10,11  > {output} 2> {log}"
           
    rule get_aligned_reads:
        input:
            aligned_reads="results/QC/{sample}-{unit}.aligned.txt",
            canonical_ids="resources/canonical_ids.csv",
        output:
            aligned_read_names="results/aligned_reads/{sample}-{unit}.read_names.txt",
        params:
            samples="results/kallisto_cds/{sample}-{unit}",
        log:
            "results/logs/aligned_read_files/{sample}-{unit}.log",
        conda:
            "../envs/QC.yaml",
        script:
            "../scripts/get_aligned_reads.py"
    
    rule get_canonical_fastq:
        input:
            get_all_fastq="results/trimmed/{sample}-{unit}.fastq.gz",
            aligned_read_names="results/aligned_reads/{sample}-{unit}.read_names.txt",
        output:
            canonical_fastq="results/canonical_reads/{sample}-{unit}.fastq",
        log:
            "results/logs/canonical_reads/{sample}-{unit}.canonical_reads.log",
        conda:
            "../envs/canonical_reads.yaml"
        shell:
            "seqtk subseq {input.get_all_fastq} {input.aligned_read_names} > {output.canonical_fastq}"  

    rule get_read_dist:
        input:
            aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples()),
        output:
            report(
                "results/QC/QC_plot.html",
                caption="../report/plot-QC.rst",
                category="QC",
            ),
        params:
            samples=expand("results/kallisto_cds/{unit.sample}-{unit.unit}",unit=units.itertuples()),
            read_length="results/stats/max-read-length.json",
        log:
            "results/logs/QC/QC_plot.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/get_histogram.py"

if config["experiment"]["3-prime-rna-seq"]["plot-qc"] != "all":

    rule get_ind_transcript_histograms:
        input:
            aligned_file=expand("results/QC/{unit.sample}-{unit.unit}.aligned.txt", unit=units.itertuples()),
        output:
            report(
                "results/QC/{transcripts}.QC_plot.html",
                caption="../report/plot-QC.rst",
                category="QC",
            ),
        params:
            transcripts =" ".join(config["experiment"]["3-prime-rna-seq"]["plot-qc"]),
            read_length="results/stats/max-read-length.json",
        log:
            "results/logs/QC/{transcripts}.QC_plot.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/histo_for_single_trans.py"