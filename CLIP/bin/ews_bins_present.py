#!/usr/bin/env python

import sys

acc = 0

thres = int(sys.argv[1])

for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')

    if cols[0] == 'chrom':
        continue

    if int(cols[4]) >= thres and int(cols[5]) >= thres:
        acc += 1

print '%s bins with EWS reads' % acc
