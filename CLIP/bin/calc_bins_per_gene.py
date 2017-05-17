#!/usr/bin/python
import sys
counter = {}
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip().split('\t')
    if 'count' in cols[4]:
        continue

    gene = cols[11]
    if not gene in counter:
        counter[gene] = 1
    else:
        counter[gene] += 1

for gene in counter:
    print counter[gene]
