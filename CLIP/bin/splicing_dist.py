#!/usr/bin/env python
#
# Splicing-site distance
# For each read (BED file), determine the distance to the nearest splice site (BED)
# Input is a BED file containing valid splice sites (and genes). Splice sites are mapped
# first to genes, then the reads are mapped to genes and the nearest enriched splice site 
# is found
#

import sys

tlss = {}
genes = {}

class Gene(object):
    
    def __init__(self, chrom, start, end, name, strand):
        self.splice_events = []
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
        pass

    def add_splice_event(self, start, end):
        self.splice_events.append((start, end))


binsize = 10000

if len(sys.argv) != 4:
    sys.stderr.write('Error - missing arguments!\nUsage: %s splice-event.bed genes.bed reads.bed\n\n' % sys.argv[0])
    sys.exit(1)


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

        gene = Gene(chrom, start, end, name, strand)

        for i in xrange(startbin, endbin+1) :
            if not (chrom, strand, i) in genes:
                genes[(chrom, strand, i)] = []

            genes[(chrom, strand, i)].append(gene)

#read splicing file (BED)
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

        mygene = None

        for i in xrange(startbin, endbin+1) :
            if mygene:
                break
            if (chrom, strand, i) in genes:
                for gene in genes[(chrom, strand, i)]:
                    if gene.start < start < gene.end or gene.start < end < gene.end:
                        mygene = gene
                        break

        if mygene:
            mygene.add_splice_event(start, end)


skipped = 0
# read in READ positions (BED)
with open(sys.argv[3]) as f:
    for line in f:
        if not line.strip() or (line and line[0] == '#'):
            continue

        cols = line.strip().split('\t')
        chrom, start, end, name, strand = cols[0], int(cols[1]), int(cols[2]), cols[3], cols[5]

        startbin = start / binsize
        endbin = end / binsize

        dist = None
        found_gene = None

        mygene = None

        for i in xrange(startbin, endbin+1) :
            if mygene:
                break
            if (chrom, strand, i) in genes:
                for gene in genes[(chrom, strand, i)]:
                    if gene.start < start < gene.end or gene.start < end < gene.end:
                        mygene = gene
                        break

        if mygene:
            best_dist = None
            best_up_down = None
            for spl_start, spl_end in mygene.splice_events:
                if spl_start <= start <= spl_end or spl_start <= end <= spl_end or start <= spl_start <= end or start <= spl_end <= end:
                    best_dist = 0
                    best_up_down = 'down'

                elif end < spl_start:
                    dist = spl_start - end
                    if best_dist is None or dist < best_dist:
                        best_dist = dist
                        best_up_down = 'up'
                elif start > spl_end:
                    dist = start - spl_end
                    if best_dist is None or dist < best_dist:
                        best_dist = dist
                        best_up_down = 'down'
                else:
                    print "not sure..."
                    print "splice: %s-%s" % (spl_start, spl_end)
                    print "read  : %s-%s" % (start, end)

            if best_dist is not None:
                if best_up_down == 'up':
                    print '%s\tSPLICE\t-%s' % (name, best_dist)
                else:
                    print '%s\tSPLICE\t%s' % (name, best_dist)
            else:
                print '%s\tSPLICE\t100000' % (name, )
                pass  # No spliced site in this gene

        else:
            print '%s\tSPLICE\t100001' % (name, )
            pass  # read not within a gene

