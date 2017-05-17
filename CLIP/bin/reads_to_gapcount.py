#!/usr/bin/env python
#
# Takes a SAM (samtools view) file as stdin, an enriched read list as argv[1]
#
# Write count files containing # gapped reads, # ungapped reads for both reads in and not in the enriched list
#

import sys

reads = set()
with open(sys.argv[1]) as f:
	for line in f:
		reads.add(line.strip())


enriched_total = 0
unenriched_total = 0

enriched_gapped = 0
unenriched_gapped = 0


for line in sys.stdin:
	cols = line.strip().split('\t')

	if cols[0] in reads:
		enriched_total += 1
	else:
		unenriched_total += 1

	cigar = []
	clen = ''
	for s in cols[5]:
		if s in '0123456789':
			clen += s
		else:
			cigar.append((s, int(clen)))
			clen = ''

	gapped = False
	for op, oplen in cigar:
		if op == 'N':
			gapped = True

	if gapped:
		if cols[0] in reads:
			enriched_gapped += 1
		else:
			unenriched_gapped += 1

print 'enriched-gapped\t%s' % enriched_gapped;
print 'enriched-ungapped\t%s' % (enriched_total-enriched_gapped);
print 'unenriched-gapped\t%s' % unenriched_gapped;
print 'unenriched-ungapped\t%s' % (unenriched_total-unenriched_gapped);

