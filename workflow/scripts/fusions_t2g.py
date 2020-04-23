import pickle
import pandas
import collections

# load dictionary mapping transcript id to:
# gene id, gene name, chromosome, start, and end
# transcript/gene version information included 
t2g = pickle.load(open(snakemake.input["t2g"], "rb"))
# dictionary mapping fusion to gene ids and gene names
fusions_t2g = {}
# dictionary mapping sequence to unique gene pairs
unique_gene_pairs = collections.defaultdict(set)
filtered = open(snakemake.output["filtered_fusions"], "w")

genes = ''
lastline = ''
with open(snakemake.input['fusions'], "r") as f:
    for l in f:
        if l[0]=='>':
            lastline = l
            l = l.replace("\n", "")
            t1, _, t2, _ = l[1:].split('_')
            genes = '_'.join([t2g[t][0] for t in (t1,t2)])
            names = '_'.join([t2g[t][1] for t in (t1,t2)])
            fusions_t2g[l[1:]] = (genes, names) # map transcript to genes and names
        elif len(l)>1:
            seq = l.replace("\n", "")
            if genes in unique_gene_pairs[seq]:
                # skip to next fasta entry if gene pair for sequence already recorded
                continue
            else:
                unique_gene_pairs[seq].add(genes)
                # write the first occuring gene pair fusion 
                filtered.write(lastline) # write header
                filtered.write(l) # write sequence

filtered.close()

def write_tsv(path, r):
    f = open(path, "w")
    for tid in r:
        f.write(tid+"\t"+"\t".join(str(x) for x in r[tid])+'\n')
    f.close()

# write fusion t2g as tsv
write_tsv(snakemake.output["fusions_t2g"], fusions_t2g)
# just keeping the gene id and gene name from the reference t2g
combined_t2g = {t:(v[0], v[1])for t, v in t2g.items()}
# combine transcript t2g and fusion t2g 
combined_t2g.update(fusions_t2g)
# write transcript and fusions t2g as tsb
write_tsv(snakemake.output["combined_t2g"], combined_t2g)
