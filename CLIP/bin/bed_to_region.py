#!/usr/bin/env python

import sys

for line in sys.stdin:
    cols = line.strip().split('\t')
    print '%s:%s-%s' % (cols[0], cols[1], cols[2])
