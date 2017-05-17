#!/usr/bin/env python
'''
Remove bins with no EWS1/EWS2/IgG/or PolyA reads
'''
import sys

for line in sys.stdin:
    if not line.strip():
        continue
    if line[0] == '#':
        sys.stdout.write(line)
        continue
    cols = line.strip('\n').split('\t')
    #   IgG                Control polyA       EWS2                EWS1
    if cols[-1] == '0' and cols[-2] == '0' and cols[-3] == '0' and cols[-4] == '0':
        continue
    sys.stdout.write(line)
