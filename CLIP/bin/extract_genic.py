#!/usr/bin/python
# Input file: EWS.enriched.bins.txt
import sys
counter = {}
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip().split('\t')
    if len(cols) <= 12:
        continue

    reg = cols[12]
    if not reg or reg == 'region_used':
        continue
    if not reg in counter:
        counter[reg] = 1
    else:
        counter[reg] += 1

for reg in sorted(counter):
    print '%s\t%s' % (reg, counter[reg])
