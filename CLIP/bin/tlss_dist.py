#!/usr/bin/env python
#
# Translational stop site (tlss) distance
# For each read (BED file), determine the distance to the nearest 3' UTR translational stop site
# Input is a BED file containing all 3'UTRs - distance calculated is from TLSS.
#
# Matches before the TLSS are negative, w/in 3' UTR are positive

import sys

tlss = {}
genes = {}

binsize = 10000

if len(sys.argv) != 4:
    sys.stderr.write('Error - missing arguments!\nUsage: %s tlss.bed genes.bed reads.bed\n\n' % sys.argv[0])
    sys.exit(1)

#read tlss file (BED)
with open(sys.argv[1]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue
        cols = line.strip().split('\t')

        chrom = cols[0]
        strand = cols[5]
        name = cols[3].split('/')[0]

        start = int(cols[1])
        end = int(cols[2])
        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1) :
            if not name in tlss:
                tlss[name] = []

            tlss[name].append((start, end))


#read gene file (BED)
with open(sys.argv[2]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue
        cols = line.strip().split('\t')

        chrom = cols[0]
        strand = cols[5]
        name = cols[3].split('/')[0]

        start = int(cols[1])
        end = int(cols[2])
        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1) :
            if not (chrom, strand, i) in genes:
                genes[(chrom, strand, i)] = []

            genes[(chrom, strand, i)].append((start, end, name))


skipped = 0
# read in READ positions (BED)
with open(sys.argv[3]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue

        cols = line.strip().split('\t')
        chrom, start, end, name, strand = cols[0], int(cols[1]), int(cols[2]), cols[3], cols[5]

        binspan = 0
        dist = None
        found_gene = False

        while not found_gene and binspan < 4:
            binspan += 1

            dist = None
            site = None

            #print "test %s:%s-%s" % (chrom, start, end)

            startbin = start / binsize
            endbin = end / binsize

            startbin -= binspan
            endbin += binspan

            #print "test %s:%s-%s (%s,%s)" % (chrom, start, end, startbin, endbin)

            best_tlss_dist = None
            best_upstream = False

            # sys.stderr.write("(%s,%s,%s, %s)\n" % (chrom, start, end, binspan))
            for i in xrange(startbin, endbin+1):
                if (chrom, strand, i) in genes: 
                    for gene_start, gene_end, gene_name in genes[(chrom, strand, i)]:
                        if gene_start <= start <= gene_end or gene_start <= end <= gene_end:
                            found_gene = True
                            if gene_name in tlss:
                                # sys.stderr.write(" => Gene: %s (%s, %s)\n" % (gene_name, gene_start, gene_end))
                                for utr_start, utr_end in tlss[gene_name]:
                                    # sys.stderr.write(" => UTR (%s, %s)\n" % (utr_start, utr_end))
                                    if strand == '+':
                                        if start <= utr_start <= end:
                                            dist = 0
                                        else:
                                            dist = start - utr_start 
                                            if dist < 0:
                                                dist = end - utr_start
                                    else:
                                        if start <= utr_end <= end:
                                            dist = 0
                                        else:
                                            dist = utr_end - end
                                            if dist < 0:
                                                dist = utr_end - start

                                    if best_tlss_dist is None or abs(dist) < best_tlss_dist:
                                        best_tlss_dist = abs(dist)
                                        if dist < 0:
                                            best_upstream = True
                                        else:
                                            best_upstream = False

            if best_tlss_dist is not None:
                dist = best_tlss_dist
                if best_upstream:
                    print '%s\tTLSS\t-%s' % (name, dist)
                else:
                    print '%s\tTLSS\t%s' % (name, dist)
