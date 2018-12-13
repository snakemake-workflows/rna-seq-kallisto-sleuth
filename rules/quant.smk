rule kallisto_index:
    input:
        config["ref"]["transcriptome"]
    output:
        "kallisto/transcripts.idx"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"


def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if len(input.fq) == 1:
        extra += " --single"
        extra += (" --fragment-length {unit.fragment_len_mean} "
                  "--sd {unit.fragment_len_sd}").format(
                    unit=units.loc[
                        (wildcards.sample, wildcards.unit)])
    else:
        extra += " --fusion"
    return extra


rule kallisto_quant:
    input:
        fq=get_fastqs,
        idx="kallisto/transcripts.idx"
    output:
        directory("kallisto/{sample}-{unit}")
    params:
        extra=kallisto_params
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output} "
        "{params.extra} {input.fq}"
