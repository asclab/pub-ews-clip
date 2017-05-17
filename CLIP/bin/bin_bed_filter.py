#!/usr/bin/env python
#
# Given a list of bins, extract all of the bins that fall within the BED file also given
#

import sys

regions = {}

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        k = (cols[0], cols[5])
        if not k in regions:
            regions[k] = set()
        regions[k].add((int(cols[1]), int(cols[2])))

with open(sys.argv[1]) as f:
    first = True
    for line in f:
        if first:
            sys.stdout.write(line)
            first = False
            continue

        cols = line.strip().split('\t')
        chrom = cols[0]
        strand = cols[3]
        start = int(cols[1])
        end = int(cols[2])
        k = (chrom, strand)
        if not k in regions:
            continue

        printed = False
        for s,e in regions[k]:
            if s <= start <= e or s <= end <= e:
                if not printed:
                    if 'UTR3' in cols[12] or 'JUNCTION' in cols[12]:
                        sys.stdout.write(line)
                        printed = True
                    continue

