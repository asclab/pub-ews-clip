#!/usr/bin/env python
#
# Inputs: Exon BED file (stranded, BED6)
#         enriched bin BED file (stranded, BED6)
#
# For each enriched bin, it finds the *smallest* exon that the bin is contained by
# Then, we calculate the location of the bin relative to the 5' donor and 3' acceptor
#
# If a bin is w/in 200bp of the donor or acceptor within the coding area of the gene,
# it will be treated as negative to the donor or plus to the acceptor. 200bp away will
# be a fixed at +/- 0.2.

import sys

bins = {}

binsize = 10000

with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        ref = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        strand = cols[5]

        exon = (ref, start, end, strand)

        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1):
            bin = (ref, i)
            if not bin in bins:
                bins[bin] = []
            bins[bin].append(exon)

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        ref = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        strand = cols[5]

        if '_INTRON' in cols[3] or 'JUNCTION' in cols[3]:
            continue

        matching_exons = set()

        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1):
            bin = (ref, i)
            if bin in bins:
                for exon in bins[bin]:
                    inref, instart, inend, instrand = exon
                    if instrand != strand:
                        continue
                    if instart < start < inend or instart < end < inend:
                        matching_exons.add(exon)

        min_exon = None
        for exon in matching_exons:
            if min_exon is None or (exon[2] - exon[1]) < (min_exon[2] - min_exon[1]):
                min_exon = exon

        if min_exon is None:
            continue

        # if start < min_exon[1] < end:
        #     # junction spans the donor site
        #     donor_dist = 0.0
        #     acceptor_dist = 1.0
        #     donor_absdist = 0
        #     acceptor_absdist = (min_exon[2] - min_exon[1])

        # elif start < min_exon[2] < end:
        #     # junction spans the acceptor site
        #     donor_dist = 1.0
        #     acceptor_dist = 0.0
        #     donor_absdist = (min_exon[2] - min_exon[1])
        #     acceptor_absdist = 0
        # else:

        # bin is w/in the exon
        acceptor_dist = float(start - min_exon[1]) / (min_exon[2] - min_exon[1])
        donor_dist = float(min_exon[2] - end) / (min_exon[2] - min_exon[1])

        acceptor_absdist = start - min_exon[1]
        donor_absdist = min_exon[2] - end

        if strand == '-':
            donor_dist, acceptor_dist = acceptor_dist, donor_dist
            donor_absdist, acceptor_absdist = acceptor_absdist, donor_absdist

        cols.append(str(-donor_dist))
        cols.append(str(acceptor_dist))
        cols.append(str(-donor_absdist))
        cols.append(str(acceptor_absdist))
        print '\t'.join(cols)
