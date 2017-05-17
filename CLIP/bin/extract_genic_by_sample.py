#!/usr/bin/python
# Input file: EWS.enriched.bins.txt
import sys
counter1 = {}
counter2 = {}
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip().split('\t')
    if len(cols) <= 12:
        continue

    reg = cols[12]

    if reg == 'region_used':
        continue

    if not reg:
        reg = 'INTERGENIC'

    if not reg in counter1:
        counter1[reg] = 0
        counter2[reg] = 0

    counter1[reg] += int(cols[4])
    counter2[reg] += int(cols[5])
    

print "region\tEWS1\tEWS2"
for reg in sorted(counter1):
    print '%s\t%s\t%s' % (reg, counter1[reg], counter2[reg])
