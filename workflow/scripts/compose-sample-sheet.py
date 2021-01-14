samples_ = snakemake.params.units[["sample", "unit"]].merge(snakemake.params.samples, on="sample")
samples_["sample"] = samples_.apply(
    lambda row: "{}-{}".format(row["sample"], row["unit"]), axis=1
)
samples_["path"] = kallisto_output
del samples_["unit"]
samples_.to_csv(output[0], sep="\t")