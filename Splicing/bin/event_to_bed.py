#!/usr/bin/env python
import sys

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')
    if cols[0] == 'event':
        continue

    # only dealing with events with multiple junctions
    if ';' in cols[0]:

        if len(cols) < 10 or not cols[9]:
            continue

        genome_span = cols[1]
        strand = cols[2]

        chrom, se = genome_span.split(':')
        start, end = se.split('-')

        print '%s\t%s\t%s\t%s\t%s\t%s' % (chrom, start, end, cols[9], 0, strand)

