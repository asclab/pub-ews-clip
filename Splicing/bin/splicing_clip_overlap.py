#!/usr/bin/env python
#
# Take a list of spliced genes (argv[1]) and a BED file for enriched bins (argv[2])
# and output a list for each spliced gene of: gene_name, is_clip_bound, clip_regions
#

import sys

clip_genes = {}

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        gene,region = cols[3].split(';')
        if not gene in clip_genes:
            clip_genes[gene] = set()

        for reg in region.split(','):
            clip_genes[gene].add(reg)

overlap = 0
print "Gene name\tclip binding in gene (Y/N)\tclip gene region (UTR/CODING)"
with open(sys.argv[1]) as f:
    for line in f:
        gene = line.strip()
        sys.stdout.write('%s\t' % gene)
        if gene in clip_genes:
            overlap += 1
            sys.stdout.write('Y\t')
            sys.stdout.write(','.join(clip_genes[gene]))
        else:
            sys.stdout.write('N\t')

        sys.stdout.write('\n')

