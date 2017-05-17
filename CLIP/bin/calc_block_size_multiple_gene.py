#!/usr/bin/env python

import sys

gene_blocks = {}

for line in sys.stdin:
    cols = line.strip().split('\t')
    start = int(cols[1])
    end = int(cols[2])

    gene = cols[3].split(';')[0]

    if not gene in gene_blocks:
        gene_blocks[gene] = []

    gene_blocks[gene].append(end-start)

for gene in gene_blocks:
    if len(gene_blocks[gene]) > 1 or gene_blocks[gene][0] > 50:
        for block in gene_blocks[gene]:
            print block

