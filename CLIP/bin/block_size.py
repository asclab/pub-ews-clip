#!/usr/bin/env python

import sys

for line in open(sys.argv[1]):
    cols = line.strip().split('\t')
    gene,region = cols[3].split(';')
    count = cols[4]

    print '%s\t%s\t%s' % (gene, region, count)
