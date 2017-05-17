#!/usr/bin/env python
import sys

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')
    if ';' in cols[0]:
        # this means that there is more than one junction in the event!

        for junc in cols[0].split(';'):
            chrom, se = junc.split(':')
            start, end = se.split('-')

            if start != end:
                print junc
