# An example collection of Snakemake rules imported in the main Snakefile.
rule index:
    input: FASTA
    output: INDEX
    log:
        "logs/kallisto/quant/{params.reference}/index.log"
    params:
        reference=REFERENCE,
        k=config['kallisto']['k']
    envmodules:
        config["envmodules"]["kallisto"]
    shell:
        "kallisto index -k {params.k} -i {output} {input} 2> {log}"

rule kallisto_quant:
    input:
        idx=INDEX,
        fastq=get_fastqs
    output:
        "kallisto/quant/{reference}/{sample}/abundance.h5",
        "kallisto/quant/{reference}/{sample}/abundance.tsv",
        "kallisto/quant/{reference}/{sample}/run_info.json",
        "kallisto/quant/{reference}/{sample}/fusion.txt"
    log: 
        "logs/kallisto/quant/{reference}/{sample}.log"
    resources:
        time=60, #minutes
        mem_mb=32000
    threads: config['kallisto']['threads']
    params:
        outdir="kallisto/quant/{reference}/{sample}",
        b=config['kallisto']['bootstrap']
    envmodules:
        config["envmodules"]["kallisto"]
    shell:
        "kallisto quant "
        "-i {input.idx} "
        "-o {params.outdir} "
        "--fusion "
        "-t {threads} "
        "-b {params.b} "
        "{input.fastq} 2> {log}"