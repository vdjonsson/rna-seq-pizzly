def unzip_fasta(fasta):
    if fasta.endswith('.gz'):
        return "<(zcat {})".format(fasta)
    else:
        return fasta

rule append_index:
    input:
        FASTA,
        "pizzly/merge/{reference}/{project}/output.fusions.combined.unique.filtered.fasta"
    output:
        "pizzly/merge/{reference}/{project}/transcripts_with_fusions.fasta.gz",
        "pizzly/merge/{reference}/{project}/transcripts_with_fusions.idx"
    log: 
        "logs/pizzly/merge/{reference}/{project}.log"
    resources:
        time=120, #minutes
        mem_mb=32000
    params:
        fa=unzip_fasta(FASTA),
        k=config['kallisto']['k']
    envmodules:
        config["envmodules"]["kallisto"]
    shell:
        "cat {params.fa} {input[1]} | gzip - > {output[0]} && "
        "kallisto index -k {params.k} -i {output[1]} {output[0]} 2> {log}"
    
rule kallisto_requant:
    input:
        idx="pizzly/merge/{reference}/{project}/transcripts_with_fusions.idx",
        fastq=get_fastqs
    output:
        "kallisto/requant/{reference}/{sample}/{project}/abundance.h5",
        "kallisto/requant/{reference}/{sample}/{project}/abundance.tsv",
        "kallisto/requant/{reference}/{sample}/{project}/run_info.json"
    log: 
        "logs/kallisto/requant/{reference}/{project}/{sample}.log"
    resources:
        time=240, #minutes
        mem_mb=32000
    threads: config['kallisto']['threads']
    params:
        outdir="kallisto/requant/{reference}/{sample}/{project}",
        b=config['kallisto']['bootstrap']
    envmodules:
        config["envmodules"]["kallisto"]
    shell:
        "kallisto quant "
        "-i {input.idx} "
        "-o {params.outdir} "
        "-t {threads} "
        "-b {params.b} "
        "{input.fastq} 2> {log}"