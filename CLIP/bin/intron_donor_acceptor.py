#!/usr/bin/env python
#
# Inputs: intron BED file (stranded, BED6)
#         enriched bin BED file (stranded, BED6)
#
# For each enriched bin, it finds the *smallest* intron that the bin is contained by
# Then, we calculate the location of the bin relative to the 5' donor and 3' acceptor
# as a relation to the size of the intron.
#
# For example, if a 50bp bin is 100bp away from the donor, and the intron was 1000bp long,
# it would be 0.1 away from the donor (100/1000) and -0.85 away from the acceptor
# (1000-(100+50))/1000
#
# [][][][][]--------------------------------------------------[][][][][][]
#                 XXXX
#          (0.1)                                       (-0.85)
#
# TODO:
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

        intron = (ref, start, end, strand)

        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1):
            bin = (ref, i)
            if not bin in bins:
                bins[bin] = []
            bins[bin].append(intron)

with open(sys.argv[2]) as f:
    for line in f:
        cols = line.strip().split('\t')
        ref = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        strand = cols[5]

        if not '_INTRON' in cols[3] and not 'JUNCTION' in cols[3]:
            continue

        matching_introns = set()

        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1):
            bin = (ref, i)
            if bin in bins:
                for intron in bins[bin]:
                    inref, instart, inend, instrand = intron
                    if instrand != strand:
                        continue
                    if instart < start < inend or instart < end < inend:
                        matching_introns.add(intron)

        min_intron = None
        for intron in matching_introns:
            if min_intron is None or (intron[2] - intron[1]) < (min_intron[2] - min_intron[1]):
                min_intron = intron

        if min_intron is None:
            continue

        if start < min_intron[1] < end:
            # junction spans the donor site
            donor_dist = 0.0
            acceptor_dist = 1.0
            donor_absdist = 0
            acceptor_absdist = (min_intron[2] - min_intron[1])

        elif start < min_intron[2] < end:
            # junction spans the acceptor site
            donor_dist = 1.0
            acceptor_dist = 0.0
            donor_absdist = (min_intron[2] - min_intron[1])
            acceptor_absdist = 0
        else:
            # bin is w/in the intron
            donor_dist = float(start - min_intron[1]) / (min_intron[2] - min_intron[1])
            acceptor_dist = float(min_intron[2] - end) / (min_intron[2] - min_intron[1])

            donor_absdist = start - min_intron[1]
            acceptor_absdist = min_intron[2] - end

        if strand == '-':
            donor_dist, acceptor_dist = acceptor_dist, donor_dist
            donor_absdist, acceptor_absdist = acceptor_absdist, donor_absdist

        cols.append(str(donor_dist))
        cols.append(str(-acceptor_dist))
        cols.append(str(donor_absdist))
        cols.append(str(-acceptor_absdist))
        print '\t'.join(cols)
