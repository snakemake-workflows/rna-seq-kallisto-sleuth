
kallisto_output = expand("kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples())


#shell.prefix("LD_LIBRARY_PATH=/home/johannes/scms/snakemake-workflows/rna-seq-kallisto-sleuth/.test/.snakemake/conda/22006fc4/lib/R/library/Rhdf5lib/lib/")


rule compose_sample_sheet:
    input:
        kallisto_output
    output:
        "sleuth/samples.tsv"
    group: "sleuth"
    run:
        samples_ = units[["sample", "unit"]].merge(samples)
        samples_["sample"] = samples_.apply(lambda row: "{}-{}".format(row["sample"], row["unit"]), axis=1)
        samples_["path"] = kallisto_output
        del samples_["unit"]
        samples_.to_csv(output[0], sep="\t")


rule sleuth_init:
    input:
        kallisto=kallisto_output,
        samples="sleuth/samples.tsv"
    output:
        "sleuth/all.rds"
    conda:
        "../envs/sleuth.yaml"
    group: "sleuth"
    script:
        "../scripts/sleuth-init.R"
