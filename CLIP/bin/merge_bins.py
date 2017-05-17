#!/usr/bin/env python
#
# Combine adjacent bins
#
# Input is the filterd bin list
# Output is a BED file
#
#
#chrom  | start       | end         | strand   | EWS1.count   | EWS2.count   | HeLa.count   | IgG.count   | gene_id   | genic_region     | gene_name   | gene_used   | region_used      | ews1_gene_total    | ews2_gene_total    | control_gene_total     | control_present_total
#chr1   | 1327900     | 1327950     | -        | 10           | 12           | 227          | 0           | 81669     | UTR3             | CCNL2       | CCNL2       | UTR3             | 315                | 398                | 186320                 | 44668                


import sys

plus_chrom = ''
plus_start = 0
plus_end = 0
plus_names = set()
plus_regions = set()
plus_count = 0

minus_chrom = ''
minus_start = 0
minus_end = 0
minus_names = set()
minus_regions = set()
minus_count = 0

for line in sys.stdin:
    if line[0] == '#':
        sys.stdout.write(line)

    cols = line.strip().split('\t')
    strand = cols[3]

    if strand not in '+-':
        continue

    chrom = cols[0]
    start = int(cols[1])
    end = int(cols[2])
    gene = cols[11]
    region = cols[12]

    if strand == '+':
        if plus_chrom == chrom and plus_end == start:
            # extend...
            plus_end = end
            plus_names.add(gene)
            plus_regions.add(region)
            plus_count += 1
        else:
            if plus_chrom:
                out = [plus_chrom, str(plus_start), str(plus_end), '%s;%s' % (','.join(plus_names), ','.join(plus_regions)), str(plus_count), '+']
                sys.stdout.write('%s\n' % '\t'.join(out))
            
            plus_chrom = chrom
            plus_start = start
            plus_end = end
            plus_names = set()
            plus_regions = set()
            plus_names.add(gene)
            plus_regions.add(region)
            plus_count = 1
    elif strand == '-':
        if minus_chrom == chrom and minus_end == start:
            # extend...
            minus_end = end
            minus_names.add(gene)
            minus_regions.add(region)
            minus_count += 1
        else:
            if minus_chrom:
                out = [minus_chrom, str(minus_start), str(minus_end), '%s;%s' % (','.join(minus_names), ','.join(minus_regions)), str(minus_count), '-']
                sys.stdout.write('%s\n' % '\t'.join(out))

            minus_chrom = chrom
            minus_start = start
            minus_end = end
            minus_names = set()
            minus_regions = set()
            minus_names.add(gene)
            minus_regions.add(region)
            minus_count = 1

out = [plus_chrom, str(plus_start), str(plus_end), '%s;%s' % (','.join(plus_names), ','.join(plus_regions)), str(plus_count), '+']
sys.stdout.write('%s\n' % '\t'.join(out))
out = [minus_chrom, str(minus_start), str(minus_end), '%s;%s' % (','.join(minus_names), ','.join(minus_regions)), str(minus_count), '-']
sys.stdout.write('%s\n' % '\t'.join(out))

