#!/usr/bin/env python
#
# Take an input BED file with GENE;REGION name, and
# split it into N smaller files based on the REGION.
#
# arg1 - input file
# arg2 = template output name

import sys

outs = {}
templ = sys.argv[2].split('%')

with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        region = cols[3].split(';')[1]
        if 'JUNCTION' in region:
            region = 'JUNCTION'
        
        if not region in outs:
            outs[region] = open('%s%s%s' % (templ[0], region, templ[1]), 'w')

        outs[region].write(line)
            
for region in outs:
    outs[region].close()
