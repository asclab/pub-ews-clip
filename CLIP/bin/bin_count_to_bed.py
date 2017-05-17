#!/usr/bin/env python

import sys
i=0
for line in sys.stdin:
    if line[0] == '#':
        continue
    
    i += 1

    cols = line.strip('\n').split('\t')
    if cols[4] != '0':
        sys.stdout.write('%s\n' % '\t'.join([cols[0], cols[1], cols[2], 'bed%s' % i, cols[4], cols[3]]))
