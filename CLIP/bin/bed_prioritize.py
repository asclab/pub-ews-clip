#!/usr/bin/env python
#
# Input is a BED file with genes and read counts as the score.
# Take in the full list, write out at least one entry for each gene,
# prioritized by the highest read-count/score. Continue until N number
# of lines have been written.
#

import sys

genes = {}

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        gene = cols[3].split(';')[0]
        score = int(cols[4])
    
        if not gene in genes:
            genes[gene] = []
        genes[gene].append((score, line))
        genes[gene].sort(reverse=True)

written_count = 0
written_genes = set()

while written_count < int(sys.argv[1]):
    best_score = 0
    best_gene = ''

    for gene in genes:
        if genes[gene]:
            if not gene in written_genes:
                sys.stdout.write(genes[gene][0][1])
                genes[gene] = genes[gene][1:]
                written_genes.add(gene)
                written_count += 1
            else:
                if genes[gene][0][0] > best_score:
                    best_score = genes[gene][0][0]
                    best_gene = gene

    if best_score:
        sys.stdout.write(genes[best_gene][0][1])
        genes[best_gene]=genes[best_gene][1:]
        written_count += 1
        written_genes.add(best_gene)

        
