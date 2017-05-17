#!/usr/bin/env python

import sys

min_count = int(sys.argv[1])

region_buf = {}
gene_counts = {}
gene_names = {}

print "## bins counts by region for genes with > %s reads" % min_count
for line in sys.stdin:
    if line[0] == '#':
        continue

    cols = line.strip('\n').split('\t')

    if not cols[5] or cols[4] == 'count':
        continue

    count = int(cols[4])
    regions = cols[6].split(',')
    geneids = cols[5].split(',')
    genenames = cols[7].split(',')

    for geneid, region, genename in zip(geneids, regions, genenames):
        if '_ANTI' in region or region == 'INTERGENIC':
            continue

        if not geneid in gene_counts:
            gene_counts[geneid] = 0
            region_buf[geneid] = {}
            gene_names[geneid] = genename

        gene_counts[geneid] += count

        if not region in region_buf[geneid]:
            region_buf[geneid][region] = 0

        region_buf[geneid][region] += 1


totals = {}

for geneid in gene_counts:
    if gene_counts[geneid] > min_count:
        for region in region_buf[geneid]:
            if not region in totals:
                totals[region] = 0
            totals[region] += region_buf[geneid][region]

for region in totals:
    print '%s\t%s' % (region, totals[region])

print ""

for geneid in gene_counts:
    if gene_counts[geneid] > min_count:
        print gene_names[geneid]
