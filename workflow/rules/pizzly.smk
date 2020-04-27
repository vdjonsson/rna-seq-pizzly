localrules: flatten, parse_gtf

rule pizzly:
    input:
        "kallisto/quant/{reference}/{sample}/fusion.txt"
    output:
        "pizzly/fusion/{reference}/{sample}/output.json",
        "pizzly/fusion/{reference}/{sample}/output.fusions.fasta"
    log: 
        "logs/pizzly/fusion/{reference}/{sample}.log"
    params:
        k=config['kallisto']['k'],
        insert_size=config['pizzly']['insert_size'],
        align_score=config['pizzly']['align_score'],
        outdir="results/pizzly/{reference}/{sample}"
    resources:
        time=120, #minutes
        mem_mb=32000
    envmodules:
        config["envmodules"]["pizzly"]
    shell:
        "pizzly "
        "-k {params.k} "
        "--gtf {GTF} "
        "--align-score {params.align_score} "
        "--insert-size {params.insert_size} "
        "--fasta {FASTA} "
        "--output {params.outdir}/output "
        "{input} 2> {log}"
    
rule flatten:
    input:
        "pizzly/fusion/{reference}/{sample}/output.json"
    output:
        "pizzly/fusion/{reference}/{sample}/genetable.tsv"
    script:
        # from https://github.com/pmelsted/pizzly
        "../scripts/flatten.py"
