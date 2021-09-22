rule init_isoform_switch:
    input:
        kallisto=kallisto_output,
        designmatrix="results/sleuth/{model}.designmatrix.rds",
        fasta="resources/transcriptome.cdna.fasta",
        gtf="resources/genome.gtf",
    output:
        rds="results/diffsplice/{model}.diffsplice.rds",
        seqs=expand(
            "results/diffsplice/{{model}}.sequences/isoformSwitchAnalyzeR_isoform_{type}.fasta",
            type=["AA", "nt"],
        ),
    log:
        "logs/diffsplice/{model}.init.log",
    conda:
        "../envs/isoform-switch-analyzer.yaml"
    params:
        model=get_model,
        samples=[
            "{unit.sample}-{unit.unit}".format(unit=unit)
            for unit in units.itertuples()
        ],
        seq_dir=lambda _, output: os.path.dirname(output.seqs[0]),
        min_effect_size=config["diffsplice"]["min_effect_size"],
        fdr=config["diffsplice"]["fdr"],
    script:
        "../scripts/isoform-switch-analysis-init.R"


rule calculate_protein_domains:
    input:
        fasta="results/diffsplice/{model}.sequences/isoformSwitchAnalyzeR_isoform_AA.fasta",
        pfam=rules.convert_pfam.output,
        pfam_dat="resources/pfam/Pfam-A.hmm.dat",
        pfam_hmm="resources/pfam/Pfam-A.hmm",
    output:
        "results/diffsplice/{model}.pfam",
    params:
        pfam_dir=lambda wildcards, input: os.path.dirname(input.pfam[0]),
    log:
        "logs/diffsplice/{model}.pfam.log",
    conda:
        "../envs/pfam.yaml"
    threads: 2
    shell:
        "pfam_scan.pl -fasta {input.fasta} -dir {params.pfam_dir} > {output} 2> {log}"


rule calculate_coding_potential:
    input:
        fasta="results/diffsplice/{model}.sequences/isoformSwitchAnalyzeR_isoform_nt.fasta",
        cpat_model="resources/cpat.logit.RData",
        hexamers="resources/cpat.hexamers.tsv",
    output:
        "results/diffsplice/{model}.cpat.tsv",
    conda:
        "../envs/cpat.yaml"
    log:
        "logs/diffsplice/{model}.cpat.log",
    shell:
        "cpat.py -g {input.fasta} -d {input.cpat_model} -x {input.hexamers} -o {output} 2> {log}"


rule annotate_isoform_switch:
    input:
        # signalP, iuPred2A, and netsufp-2 are not distributable via Conda
        rds="results/diffsplice/{model}.diffsplice.rds",
        pfam="results/diffsplice/{model}.pfam",
        cpat="results/diffsplice/{model}.cpat.tsv",
    output:
        plots_with=report(
            directory("results/plots/diffsplice/{model}/with_consequences"),
            category="Differential splicing  plots (with consequences)",
            subcategory="{model}",
            patterns=["{no}_switch_plot_{geneid}_aka_{gene}.pdf"],
            caption="../report/diffsplice-plot.rst",
        ),
        plots_without=report(
            directory("results/plots/diffsplice/{model}/without_consequences"),
            category="Differential splicing plots (without consequences)",
            subcategory="{model}",
            patterns=["{no}_switch_plot_{geneid}_aka_{gene}.pdf"],
            caption="../report/diffsplice-plot.rst",
        ),
        table=report(
            "results/tables/diffsplice/{model}.diffsplice.tsv",
            category="Differential splicing",
            subcategory="{model}",
            caption="../report/diffsplice.rst",
        ),
    log:
        "logs/diffsplice/{model}.annotate.log",
    params:
        coding_cutoff=config["diffsplice"]["coding_cutoff"],
        remove_noncoding_orfs=config["diffsplice"]["remove_noncoding_orfs"],
        plotdir=lambda _, output: os.path.dirname(output.plots_with),
        min_effect_size=config["diffsplice"]["min_effect_size"],
        fdr=config["diffsplice"]["fdr"],
    conda:
        "../envs/isoform-switch-analyzer.yaml"
    script:
        "../scripts/isoform-switch-analysis-annotate.R"
