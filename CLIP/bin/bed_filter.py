#!/usr/bin/env python

import sys

genes = set()

with open(sys.argv[1]) as f:
    for line in f:
        gene = line.strip()
        genes.add(gene)

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        if cols[3].split(';')[0] in genes:
            sys.stdout.write(line)

