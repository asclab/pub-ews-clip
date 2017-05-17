#!/usr/bin/env python
import sys

genes = set()

with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        if cols[11] and cols[11] != 'gene_used':
            genes.add(cols[11])

for gene in sorted(genes):
    print gene
