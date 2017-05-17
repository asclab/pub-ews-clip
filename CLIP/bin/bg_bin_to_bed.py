#!/usr/bin/env python

import sys


def run(fname):
    header = True
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue
            if not line.strip():
                continue
            if header:
                header = False
                continue

            cols = line.strip('\n').split('\t')

            hela_count = int(cols[6])
            if hela_count == 0:
                continue

            region = cols[12]
            if not region:
                continue

            sys.stdout.write('%s\n' % '\t'.join([cols[0], cols[1], cols[2], '%s;%s' % (cols[11], cols[12]), str(hela_count), cols[3]]))

if __name__ == '__main__':
    run(sys.argv[1])
