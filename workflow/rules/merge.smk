localrules: cat_fasta, unique_fasta, count_fasta, fusions_t2g

rule cat_fasta:
    input:
        expand("pizzly/fusion/{reference}/{sample}/output.fusions.fasta", 
            reference=REFERENCE, sample=samples['sample'])
    output:
        metadata="pizzly/merge/{reference}/{project}/samples.txt",
        fa="pizzly/merge/{reference}/{project}/output.fusions.combined.fasta"
    resources:
        time=120,
        mem_mb=8000
    shell:
        """
        cat {input} > {output.fa} 

        echo "Samples:" > {output.metadata} 
        for file in {input}
        do 
            echo $file >> {output.metadata} 
        done
        """
    
rule unique_fasta:
    input:
        "pizzly/merge/{reference}/{project}/{prefix}.fasta"
    output:
        "pizzly/merge/{reference}/{project}/{prefix}.unique.fasta"
    resources:
        time=120, #minutes
        mem_mb=8000
    shell:
        "awk '!a[$0]++' RS='>' ORS='>' {input} | awk '!/^>?$/' > {output} "

rule count_fasta:
    input:
        raw="pizzly/merge/{reference}/{project}/output.fusions.combined.fasta",
        unique="pizzly/merge/{reference}/{project}/output.fusions.combined.unique.fasta",
        filtered="pizzly/merge/{reference}/{project}/output.fusions.combined.unique.filtered.fasta"
    output:
        "pizzly/merge/{reference}/{project}/summary.txt"
    shell:
        """        
        echo "Raw  {input.raw}:" > {output} 
        {{ printf "Total number of headers: " ; grep '>' {input.raw} | wc -l ;}} >> {output} 
        {{ printf "Total number of sequences: " ; grep -v '>' {input.raw} | wc -l ;}} >> {output}
        {{ printf "Number of unique headers: " ; grep '>' {input.raw} | sort | uniq | wc -l ;}} >> {output} 
        {{ printf "Number of unique sequences: " ; grep -v '>' {input.raw} | sort | uniq | wc -l ;}} >> {output}  
        echo '-----------------------' >> {output} 
        echo "Unique {input.unique}:"  >> {output} 
        {{ printf "Total number of headers: " ; grep '>' {input.unique} | wc -l ;}} >> {output} 
        {{ printf "Total number of sequences: " ; grep -v '>' {input.unique} | wc -l ;}} >> {output} 
        {{ printf "Number of unique headers: " ; grep '>' {input.unique} | sort | uniq | wc -l ;}} >> {output} 
        {{ printf "Number of unique sequences: " ; grep -v '>' {input.unique} | sort | uniq | wc -l ;}} >> {output} 
        echo '-----------------------' >> {output} 
        echo "Filtered {input.filtered}:"  >> {output} 
        {{ printf "Total number of headers: " ; grep '>' {input.filtered} | wc -l ;}} >> {output} 
        {{ printf "Total number of sequences: " ; grep -v '>' {input.filtered} | wc -l ;}} >> {output} 
        {{ printf "Number of unique headers: " ; grep '>' {input.filtered} | sort | uniq | wc -l ;}} >> {output} 
        {{ printf "Number of unique sequences: " ; grep -v '>' {input.filtered} | sort | uniq | wc -l ;}} >> {output} 
        """

rule parse_gtf:
    input:
        GTF
    output:
        T2G
    # add pickle and pandas requirement
    script:
        "../scripts/parse_gtf.py"
        
rule fusions_t2g:
    input:
        t2g=T2G,
        fusions="pizzly/merge/{reference}/{project}/output.fusions.combined.unique.fasta"
    output:
        fusions_t2g="pizzly/merge/{reference}/{project}/fusions_t2g.tsv",
        combined_t2g="pizzly/merge/{reference}/{project}/combined_t2g.tsv",
        filtered_fusions="pizzly/merge/{reference}/{project}/output.fusions.combined.unique.filtered.fasta"
    # add pickle and pandas requirement
    script:
        "../scripts/fusions_t2g.py"

