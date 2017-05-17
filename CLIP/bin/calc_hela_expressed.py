#!/usr/bin/env python
#
# Tallies the number of HeLa reads and bins by genic region
#

import sys

region_bins = {}
region_reads = {}

priority = ['JUNCTION', 'CODING', 'UTR5', 'UTR3', 'CODING_INTRON', 'UTR5_INTRON', 'UTR3_INTRON', ]

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')
    if cols[0] == 'chrom':
        continue

    hela_reads = int(cols[6])

    if hela_reads > 0:
        regions = set(cols[9].split(','))
        region = None

        for prior in priority:
            if prior in regions:
                region = prior
                break

        if not region:
            # if we aren't in a region above (by default means we are w/in a gene and not anti-sense)
            continue

        if not region in region_bins:
            region_bins[region] = 1
        else:
            region_bins[region] += 1

        if not region in region_reads:
            region_reads[region] = hela_reads
        else:
            region_reads[region] += hela_reads

for region in sorted(region_bins):
    print '%s\t%s\t%s' % (region, region_bins[region], region_reads[region])
