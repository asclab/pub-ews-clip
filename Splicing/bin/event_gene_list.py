#!/usr/bin/env python
import sys

single = 0
novel = 0
total = 0
gene_counts = {}

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')
    if cols[0] == 'event':
        continue

    # only dealing with events with multiple junctions
    if ';' in cols[0]:
        total += 1

        if len(cols) < 10 or not cols[9]:
            novel += 1
            continue

        for gene in cols[9].split(','):
            if not gene in gene_counts:
                gene_counts[gene] = 1
            else:
                gene_counts[gene] += 1
    else:
        single += 1

sys.stderr.write('Total events (multiple-junctions): %s\n' % total)
sys.stderr.write('Single-junction events: %s\n' % single)
sys.stderr.write('Total genes: %s\n' % len(gene_counts))
sys.stderr.write('Novel genes: %s\n' % novel)

for gene in gene_counts:
    print gene
