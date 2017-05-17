#!/usr/bin/env python

#
# Given a comma-delimited list of genic regions, collapse them into a prioritized list
#

import sys

priority = ['JUNCTION', 'CODING', 'UTR5', 'UTR3', 'CODING_INTRON', 'UTR5_INTRON', 'UTR3_INTRON', ]

for line in sys.stdin:
    regions = set(line.strip().split(','))

    for prior in priority:
        if prior in regions:
            print prior
            break
