#!/usr/bin/env python

# Adapted from https://github.com/bustools/getting_started/releases/ 
# t2g.py

import sys, argparse
import pickle

def load_gtf(input, use_name=True, use_version=True):
    r = {}
    for line in input:
        if len(line) == 0 or line[0] == '#':
            continue
        l = line.strip().split('\t')
        if l[2] == 'transcript':
            info = l[8]
            d = {}
            for x in info.split('; '):
                x = x.strip()
                p = x.find(' ')
                if p == -1:
                    continue
                k = x[:p]
                p = x.find('"',p)
                p2 = x.find('"',p+1)
                v = x[p+1:p2]
                d[k] = v

            if 'transcript_id' not in d or 'gene_id' not in d:
                continue

            tid = d['transcript_id'].split(".")[0]
            gid = d['gene_id'].split(".")[0]
            if use_version:
                if 'transcript_version' not in d or 'gene_version' not in d:
                    continue
                tid += '.' + d['transcript_version']
                gid += '.' + d['gene_version']
            gname = None
            if use_name:
                if 'gene_name' not in d:
                    continue
                gname = d['gene_name']
            # position 
            chrom = l[0]
            start = l[3]
            end = l[4]
            if tid in r:
                continue

            r[tid] = (gid, gname, chrom, start, end)
    return r

def write_pickle(path, r):
    f = open(path, "wb")
    pickle.dump(r, f)
    f.close()

def write_tsv(path, r):
    f = open(path, "w")
    for tid in r:
        f.write(tid+"\t"+"\t".join(str(x) for x in r[tid])+'\n')
    f.close()

def load_pickle(path):
    return pickle.load(open(path, "rb"))

if __name__ == "__main__":
    gtf = open(snakemake.input[0], "r")
    r = load_gtf(gtf, use_name=True, use_version=True)
    write_pickle(snakemake.output[0], r)
    gtf.close()
