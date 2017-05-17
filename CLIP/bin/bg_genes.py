#!/usr/bin/env python

import sys

genes = set()

with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        if cols[0] == 'chrom':
            continue

        ews_counts = int(cols[4]) + int(cols[5])

        if ews_counts > 0 and cols[11]:
            genes.add(cols[11])

for gene in genes:
    print gene
