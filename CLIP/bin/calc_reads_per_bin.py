#!/usr/bin/python
# Input file: EWS.enriched.bins.txt
import sys
counter = {}
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip().split('\t')
    if 'count' in cols[4]:
        continue
    # just print EWS1 + EWS2
    print int(cols[4]) + int(cols[5])
