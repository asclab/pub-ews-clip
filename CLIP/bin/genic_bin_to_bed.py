#!/usr/bin/env python
'''
convert genic-annotated bins to a bed file with colors...
'''

import sys
import gzip
import itertools

color_vals = list(itertools.combinations_with_replacement(['0', '80', '128', '192', '255'], 3))
color_idx = 0
sys.stderr.write('%s\n' % ';'.join(['(%s,%s,%s)' % x for x in color_vals]));

color_vals = color_vals[1:-1]

colors = {}

first = True

genic_region_col = 9
gene_name_col = 10

sys.stdout.write('track name="bins" itemRgb="On"\n')

with gzip.open(sys.argv[1]) as f:
	for line in f:
		if line[0] == '#':
			continue

		if first:
			first = False
			continue

		cols = line.strip('\n').split('\t')
		ref = cols[0]
		start = cols[1]
		end = cols[2]
		strand = cols[3]
		name = cols[gene_name_col]
		genic_region = cols[genic_region_col]

		if not genic_region or ',' in genic_region:
			color = '255,255,255'
		else:
			if not genic_region in colors:
				colors[genic_region] = '%s,%s,%s' % (color_vals[color_idx])
				sys.stderr.write('%s => %s\n' % (genic_region, colors[genic_region]))
				color_idx += 1

			color = colors[genic_region]

		sys.stdout.write('%s\n' % '\t'.join([ref, start, end, '%s-%s' % (genic_region, name), '0', strand, start, end, color]))
