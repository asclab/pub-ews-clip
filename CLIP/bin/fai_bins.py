#!/usr/bin/env python
import sys

binsize = 50

# FAI file
with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        if '_' in cols[0]:
            continue
        num = 0
        pos = 0
        while pos < int(cols[1]):
            out = [cols[0], pos, pos+binsize, '%s_%s' % (cols[0], num), '0', '+']
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in out]))
            out = [cols[0], pos, pos+binsize, '%s_%s' % (cols[0], num), '0', '-']
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in out]))
            pos += binsize
            num += 1
