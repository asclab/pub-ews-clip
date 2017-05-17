#!/usr/bin/env python
import sys

genes = set()

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')

    if cols[0] == 'junction':
        continue

    if cols[-1]:
    	for gene in cols[-1].split(','):
        	genes.add(gene)

for gene in genes:
    print gene
