#!/usr/bin/env python

import sys


def run(f):
    header = True
    for line in f:
        if line[0] == '#':
            continue
        if not line.strip():
            continue
        if header:
            header = False
            continue

        cols = line.strip('\n').split('\t')

        sys.stdout.write('%s\n' % '\t'.join([cols[0], cols[1], cols[2], '%s;%s' % (cols[11], cols[12]), str(int(cols[4]) + int(cols[5])), cols[3]]))

if __name__ == '__main__':
    if sys.argv[1] == '-':
        run(sys.stdin)
    else:
        with open(sys.argv[1]) as f:
            run(f)
