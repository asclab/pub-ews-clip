#!/usr/bin/python
# Input file: EWS.enriched.bins.txt
import sys
ews1_count = {}
ews2_count = {}
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip().split('\t')

    if len(cols) <= 12:
        continue

    if 'count' in cols[4]:
        continue

    ews1 = int(cols[4])
    ews2 = int(cols[5])
    reg = cols[12]

    if not reg or reg == 'region_used':
        continue

    if not reg in ews1_count:
        ews1_count[reg] = 0

    if not reg in ews2_count:
        ews2_count[reg] = 0

    ews1_count[reg] += ews1
    ews2_count[reg] += ews2

print "Region\tEWS-1\tEWS-2"
for reg in sorted(ews1_count):
    print '%s\t%s\t%s' % (reg, ews1_count[reg], ews2_count[reg])
