
kallisto_output = expand("results/requant/{reference}/{sample}/{project}", 
    reference=REFERENCE, sample=samples['sample'], project=PROJECT)

localrules: compose_sample_sheet

rule compose_sample_sheet:
    input:
        kallisto_output
    output:
        "sleuth/{reference}/{project}/samplesheet.tsv"
    run:
        samples_ = samples.copy()
        samples_['path'] =  kallisto_output
        samples_.to_csv(output[0], sep="\t", index=False)

def get_model(wildcards):
    if wildcards.model == "all":
        return {"full": None}
    return config["sleuth"]["models"][wildcards.model]

rule check_h5:
    input: 
        expand("kallisto/requant/{reference}/{sample}/{project}/abundance.h5", 
            reference=REFERENCE, sample=samples['sample'], project=PROJECT),
        samples="sleuth/{reference}/{project}/samplesheet.tsv"
    output:
        "sleuth/{reference}/{project}/h5_status.txt"
    conda:
        "../envs/sleuth.yaml"
    envmodules:
        config["envmodules"]['sleuth']
    resources:
        time=120, #minutes
        mem_mb=32000
    script:
        "../scripts/check_h5.R"

