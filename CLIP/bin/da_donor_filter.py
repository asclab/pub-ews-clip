#!/usr/bin/env python

import sys

bed = False

if sys.argv[1] == '-bed':
    bed = True
    fname = sys.argv[2]
else:
    fname = sys.argv[1]

with open(fname) as f:
    for line in f:
        cols = line.strip().split('\t')
        donor = int(cols[-2])
        acceptor = int(cols[-1])

        name = cols[3]
        region = name.split(';')[1]

        if 'UTR3' in region:
            continue
        if 'UTR5' in region:
            continue


        if abs(donor) < abs(acceptor) and abs(donor) <= 500:
            if bed:
                sys.stdout.write('%s\n' % '\t'.join(cols[:5]))
            else:
                sys.stdout.write(line)
