import sys
import json
from collections import OrderedDict

def loadJSON(fn):
    with open(fn) as f:
        JJ = json.load(f,object_pairs_hook=OrderedDict)
    return JJ['genes']

def outputGeneTable(fusions, outf, filters = None):
    outf.write('\t'.join("geneA.name geneA.id geneB.name geneB.id paircount splitcount transcripts.list".split()))
    outf.write('\n')
    for gf in fusions:
        gAname = gf['geneA']['name']
        gAid   = gf['geneA']['id']
        gBname = gf['geneB']['name']
        gBid   = gf['geneB']['id']
        pairs  = str(gf['paircount'])
        split  = str(gf['splitcount'])
        txp = [tp['fasta_record'] for tp in gf['transcripts']]

        outf.write('\t'.join([gAname, gAid, gBname, gBid, pairs, split, ';'.join(txp)]))
        outf.write('\n')

outf = open(snakemake.output[0],'w')
fusions = loadJSON(snakemake.input[0])
outputGeneTable(fusions, outf)
outf.close()