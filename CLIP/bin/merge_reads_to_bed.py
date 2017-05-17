#!/usr/bin/env python
#
# Take a BED file of reads, and a distance file showing how close a read 
# is to 3'UTR, DONOR or ACCEPTOR.
#
# Output a BED file with reads that belong to one of those classes (argv[1] in ['utr3', 'acceptor', 'donor'])
#
##


import sys

if len(sys.argv) != 4 or sys.argv[1] not in ['utr3', 'donor', 'acceptor']:
    sys.stderr.write("Error! Usage: merge_reads_to_bed.py [utr3|acceptor|donor] reads.bed reads.dist.txt\n\n")
    sys.exit(1)

reads = set()

dist_exon = 500
dist_intron = 2000

with open(sys.argv[3]) as f:
    for line in f:
        name, clazz, dist = line.strip().split('\t')

        if sys.argv[1] == 'utr3' and clazz == 'UTR3':
            reads.add(name)
        if sys.argv[1] == 'donor' and clazz == 'DONOR':
            if -dist_exon <= dist <= dist_intron:
                reads.add(name)
        if sys.argv[1] == 'acceptor' and clazz == 'ACCEPTOR':
            if -dist_intron <= dist <= dist_exon:
                reads.add(name)

with open(sys.argv[2]) as f:
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        if cols[3] in reads:
            sys.stdout.write(line)
