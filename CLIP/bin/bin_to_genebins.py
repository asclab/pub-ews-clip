#!/usr/bin/env python

import sys

genes = {}

def run(fname):
    header = True
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue
            if not line.strip():
                continue
            if header:
                header=False
                continue
            

            cols = line.strip('\n').split('\t')
            if not cols[11] in genes:
                genes[cols[11]] = 1
            else:
                genes[cols[11]] += 1

    for gene in sorted(genes):
        print '%s\t%s' % (gene, genes[gene])

if __name__ == '__main__':
    run(sys.argv[1])
