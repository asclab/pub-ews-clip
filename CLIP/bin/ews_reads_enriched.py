#!/usr/bin/env python

import sys

acc = 0

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')

    if cols[0] == 'chrom':
        continue

    acc += int(cols[4]) + int(cols[5])

print '%s EWS reads in enriched bins' % acc
