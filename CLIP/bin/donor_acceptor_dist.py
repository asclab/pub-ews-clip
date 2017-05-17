#!/usr/bin/env python
#
# Donor - acceptor distance
# For each read (BED file), determine the nearest Donor or Acceptor site
#
# Donor matches in exon are negative, in the intron are positive
# Acceptor in exon are positive, intron matches are negative
#
# Test - removing reads w/in 3'UTR
#

import sys

donors = {}
acceptors = {}
tlss = {}

binsize = 10000

if len(sys.argv) != 5:
    sys.stderr.write('Error - missing arguments!\nUsage: %s donors.bed acceptors.bed tlss.bed reads.bed\n\n' % sys.argv[0])
    sys.exit(1)

#read donors file (BED)
with open(sys.argv[1]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue
        cols = line.strip().split('\t')

        chrom = cols[0]
        strand = cols[5]

        pos = int(cols[1]) if strand == '+' else int(cols[2])
        bin = pos / binsize

        if not (chrom, strand, bin) in donors:
            donors[(chrom, strand, bin)] = []

        donors[(chrom, strand, bin)].append(pos)
        #print "donors[(%s, %s, %s)].append(%s) %s:%s-%s" % (chrom, strand, bin, pos, chrom, cols[1], cols[2])


#read acceptors file (BED)
with open(sys.argv[2]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue
        cols = line.strip().split('\t')

        chrom = cols[0]
        strand = cols[5]

        pos = int(cols[2]) if strand == '+' else int(cols[1])
        bin = pos / binsize

        if not (chrom, strand, bin) in acceptors:
            acceptors[(chrom, strand, bin)] = []

        acceptors[(chrom, strand, bin)].append(pos)
        #print "acceptors[(%s, %s, %s)].append(%s) %s:%s-%s" % (chrom, strand, bin, pos, chrom, cols[1], cols[2])


### UTR3 BED file
with open(sys.argv[3]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue
        cols = line.strip().split('\t')

        chrom = cols[0]
        strand = cols[5]

        start = int(cols[1])
        end = int(cols[2])
        startbin = start / binsize
        endbin = end / binsize

        for i in xrange(startbin, endbin+1) :
            if not (chrom, strand, i) in tlss:
                tlss[(chrom, strand, i)] = []

            tlss[(chrom, strand, i)].append((start, end))


skipped = 0
# read in READ positions (BED)
with open(sys.argv[4]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue

        cols = line.strip().split('\t')
        chrom, start, end, name, strand = cols[0], int(cols[1]), int(cols[2]), cols[3], cols[5]

        binspan = 0
        dist = None

        while dist is None:
            binspan += 1

            dist = None
            site = None

            #print "test %s:%s-%s" % (chrom, start, end)

            startbin = start / binsize
            endbin = end / binsize

            startbin -= binspan
            endbin += binspan


            inUTR = False

            for i in xrange(startbin, endbin+1):
                if (chrom, strand, i) in tlss:
                    for utr_start, utr_end in tlss[(chrom, strand, i)]:
                        if utr_start <= start <= utr_end:
                            inUTR = True
                        if utr_start <= end <= utr_end:
                            inUTR = True

            if inUTR:
                print '%s\tUTR3\t0' % (name)
                dist = 0
                continue

            #print "test %s:%s-%s (%s,%s)" % (chrom, start, end, startbin, endbin)

            best_donor_dist = None
            best_acceptor_dist = None

            best_donor_upstream = False
            best_acceptor_upstream = False

            for i in xrange(startbin, endbin+1):
                if (chrom, strand, i) in donors:
                    for pos in donors[(chrom, strand, i)]:
                        if start <= pos <= end:
                            best_donor_dist = 0
                            continue
                        else:
                            if (end < pos):
                                dist = pos - end
                                if best_donor_dist is None or dist < best_donor_dist:
                                    best_donor_dist = dist
                                    best_donor_upstream = True
                            else:
                                dist = start - pos
                                if best_donor_dist is None or dist < best_donor_dist:
                                    best_donor_dist = dist
                                    best_donor_upstream = False

            for i in xrange(startbin, endbin+1):
                if (chrom, strand, i) in acceptors:
                    for pos in acceptors[(chrom, strand, i)]:
                        if start <= pos <= end:
                            best_acceptor_dist = 0
                            continue
                        else:
                            if (end < pos):
                                dist = pos - end
                                if best_acceptor_dist is None or dist < best_acceptor_dist:
                                    best_acceptor_dist = dist
                                    best_acceptor_upstream = True
                            else:
                                dist = start - pos
                                if best_acceptor_dist is None or dist < best_acceptor_dist:
                                    best_acceptor_dist = dist
                                    best_acceptor_upstream = False


            if best_donor_dist == 0 and best_acceptor_dist == 0:
                print '%s\tDONOR\t0\n%s\tACCEPTOR\t0' % (name, name)
                dist = 0
                continue

            if best_acceptor_dist is None or best_donor_dist < best_acceptor_dist:
                dist = best_donor_dist
                best_upstream = best_donor_upstream
                site = "DONOR"
            else:
                dist = best_acceptor_dist
                best_upstream = best_acceptor_upstream
                site = "ACCEPTOR"

            if dist is None:
                continue

            if strand == '+':
                if best_upstream:
                    dist = -dist
            else:
                if not best_upstream:
                    dist = -dist

            print '%s\t%s\t%s' % (name, site, dist)
