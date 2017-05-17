#!/usr/bin/env python
#
# For each bin, if the HeLa count > 0, what is the total EWS count?
# This will be used for the background of a density plot.
#

import sys

priority = ['JUNCTION', 'CODING', 'UTR5', 'UTR3', 'CODING_INTRON', 'UTR5_INTRON', 'UTR3_INTRON', ]

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')
    if cols[0] == 'chrom':
        continue

    hela_reads = int(cols[6])
    ews_reads = int(cols[4]) + int(cols[5])

    if hela_reads > 0 and cols[8]:  # must be w/in a gene...
        regions = set(cols[9].split(','))
        region = None

        for prior in priority:
            if prior in regions:
                region = prior
                break

        if not region:
            # if we aren't in a region above bail
            # (by default means we are w/in a gene and not anti-sense)
            continue

        ews_reads = int(cols[4]) + int(cols[5])
        print ews_reads
