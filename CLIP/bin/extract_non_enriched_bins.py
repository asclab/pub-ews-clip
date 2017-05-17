#!/usr/bin/env python
import sys

enriched_fname = sys.argv[1]
all_fname = sys.argv[2]

enriched = set()
enriched_genes = set()

sys.stderr.write('Loading enriched bins...\n')
with open(enriched_fname) as f:
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        if cols[0] == 'chrom':
            continue
        k = (cols[0], cols[1], cols[2], cols[3])
        enriched.add(k)
        enriched_genes.add(cols[11])
sys.stderr.write("Enriched genes: %s" % len(enriched_genes))
sys.stderr.write('Loading all bins...\n')
i=0
with open(all_fname) as f:
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        if cols[0] == 'chrom':
            sys.stdout.write('%s\n' % '\t'.join(cols[:6]))
            continue
        if cols[0] == 'chrM':
            continue
        k = (cols[0], cols[1], cols[2], cols[3])
        if not k in enriched:
            if not cols[11]:
                # not in a gene
                continue
            if cols[4] != '0' or cols[5] != '0':
                i += 1
                if i > 10000:
                    i=0
                    sys.stderr.write('%s\n' % str(k))
                sys.stdout.write('%s\n' % '\t'.join(cols))
