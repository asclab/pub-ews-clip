#!/usr/bin/env python
'''
Calculates the CLIP enrichment score for each bin using the control polyA RNAseq
counts to estimate the number of CLIP reads. This performs a Poisson test of the
actual CLIP bin counts to an estimate of:
    EWS reads in a particular region-class (coding, utrs, introns, junction, etc...):
    Expected = # EWS reads in a region-class * (# control reads in bin + 1 / # control reads in region-class )

The p-value is based upon how likley the actual count would occur if the "true
value" was the highest of the estimated rates.

Bins with a p-value higher than the threshold for either EWS1 or EWS2 will be removed
(default threshold: 0.01)
'''

import sys
import gzip
import collections

# how many reads are in each gene - gene_name is the key for these
ews1_counts = collections.defaultdict(int)
ews2_counts = collections.defaultdict(int)
control_counts = collections.defaultdict(int)

ews1_region_counts = collections.defaultdict(int)
ews2_region_counts = collections.defaultdict(int)
control_region_counts = collections.defaultdict(int)

gene_names = {}

# how many bins should belong to each gene
# gene_bin_size = {}

ews1_col = 4
ews2_col = 5
control_col = 6
gene_col = 8
gene_name_col = 10
genic_region_col = 9


def openfile(fname):
    if fname[-3:] == '.gz':
        return gzip.open(fname)
    return open(fname)


def calc_prob(fname, gene_size_fname, thres=None, min_count=5):

    # sys.stderr.write("Calculating gene bin counts\n")
    # with openfile(gene_size_fname) as f:
    #     for line in f:
    #         cols = line.strip().split('\t')
    #         start = int(cols[1])
    #         end = int(cols[2])
    #         name = cols[3]

    #         start_bin = start / 50
    #         end_bin = end / 50

    #         gene_bin_size[name] = end_bin - start_bin + 1

    header = True
    sys.stderr.write("Calculating gene totals\n")
    i = 0

    ews1_total = 0
    ews2_total = 0
    control_total = 0

    with openfile(fname) as f:
        for line in f:
            i += 1
            if not line.strip() or line[0] == '#':
                continue
            if header:
                header = False
                continue

            cols = line.strip('\n').split('\t')

            if cols[0] == 'chrM':
                continue

            if i > 100000:
                sys.stderr.write('%s:%s\n' % (cols[0], cols[1]))
                i = 0

            if not cols[gene_col]:
                continue

            genes = cols[gene_col].split(',')
            names = cols[gene_name_col].split(',')
            regions = cols[genic_region_col].split(',')

            ews1_count = int(cols[ews1_col]) if cols[ews1_col] else 0
            ews2_count = int(cols[ews2_col]) if cols[ews2_col] else 0
            control_count = int(cols[control_col]) if cols[control_col] else 0

            ews1_total += ews1_count
            ews2_total += ews2_count
            control_total += control_count

            for gene, name, region in zip(genes, names, regions):
                gene_names[gene] = name
                ews1_counts[gene] += ews1_count
                ews2_counts[gene] += ews2_count
                control_counts[gene] += control_count

                ews1_region_counts[(gene, region)] += ews1_count
                ews2_region_counts[(gene, region)] += ews2_count
                control_region_counts[(gene, region)] += control_count

    i = 0

    sys.stdout.write('# minimum EWS read count: %s\n' % min_count)
    sys.stdout.write('# EWS1 total: %s\n' % ews1_total)
    sys.stdout.write('# EWS2 total: %s\n' % ews2_total)
    sys.stdout.write('# Control (Nuclear RNA) total: %s\n' % control_total)

    ews1_factor = float(ews1_total) / control_total
    ews2_factor = float(ews2_total) / control_total

    sys.stdout.write('# HeLa-EWS1 factor: %s\n' % ews1_factor)
    sys.stdout.write('# HeLa-EWS2 factor: %s\n' % ews2_factor)

    sys.stderr.write('Calculating enrichment\n')
    header = True
    with openfile(fname) as f:
        for line in f:
            i += 1
            if not line.strip() or line[0] == '#':
                sys.stdout.write(line)
                continue

            cols = line.strip('\n').split('\t')
            if cols[0] == 'chrM':
                continue
            if i > 100000:
                sys.stderr.write('%s:%s\n' % (cols[0], cols[1]))
                i = 0

            if header:
                cols.append('gene_used')
                cols.append('region_used')

                cols.append('ews1_gene_total')
                cols.append('ews2_gene_total')
                cols.append('control_gene_total')

                cols.append('ews1_gene_region')
                cols.append('ews2_gene_region')
                cols.append('control_gene_region')

                cols.append('effective_control_region_pct')

                cols.append('ews1_expected_region')
                cols.append('ews2_expected_region')

                header = False
                sys.stdout.write('%s\n' % '\t'.join(cols))
                continue

            gene_used = ''
            region_used = ''
            ews1_gene_total = 0
            ews2_gene_total = 0
            control_gene_total = 0

            ews1_count = int(cols[ews1_col]) if cols[ews1_col] else 0
            ews2_count = int(cols[ews2_col]) if cols[ews2_col] else 0
            control_count = int(cols[control_col]) if cols[control_col] else 0

            ews1_gene_region = 0
            ews2_gene_region = 0
            control_gene_region = 0

            ews1_expected_region = 0
            ews2_expected_region = 0

            # if ews1_count < min_count and ews2_count < min_count:
            #     # require reads in  EWS1 or EWS2
            #     continue

            genes = cols[gene_col].split(',')
            genenames = cols[gene_name_col].split(',')
            regions = cols[genic_region_col].split(',')

            if genes:
                best_gene = genes[0]
                gene_used = genenames[0]
                region_used = regions[0]
                best_total = control_counts[genes[0]]
                for j, geneid in enumerate(genes):
                    if control_counts[geneid] > best_total:
                        try:
                            best_total = control_counts[geneid]
                            best_gene = geneid
                            gene_used = genenames[j] if len(genenames) > j else geneid
                            region_used = regions[j]
                        except Exception, e:
                            print e
                            print cols[0], cols[1], cols[2]
                            print genes
                            print genenames
                            print regions
                            sys.exit(1)

                ews1_gene_total = ews1_counts[best_gene]
                ews2_gene_total = ews2_counts[best_gene]
                control_gene_total = control_counts[best_gene]

                ews1_gene_region = ews1_region_counts[(best_gene, region_used)]
                ews2_gene_region = ews2_region_counts[(best_gene, region_used)]
                control_gene_region = control_region_counts[(best_gene, region_used)]

            if control_gene_region > 0:
                effective_control_region_pct = float(control_count + 1) / control_gene_region
            else:
                effective_control_region_pct = 0.0

            if ews1_count > 0:
                if ews1_gene_region > 0:
                    ews1_expected_region = effective_control_region_pct * ews1_gene_region

            if ews2_count > 0:
                if ews2_gene_region > 0:
                    ews2_expected_region = effective_control_region_pct * ews2_gene_region

            cols.append(gene_used)
            cols.append(region_used)

            cols.append(ews1_gene_total)
            cols.append(ews2_gene_total)
            cols.append(control_gene_total)

            cols.append(ews1_gene_region)
            cols.append(ews2_gene_region)
            cols.append(control_gene_region)

            cols.append(effective_control_region_pct)

            cols.append(ews1_expected_region)
            cols.append(ews2_expected_region)

            # if not thres or (ews1_pvalue < thres and ews2_pvalue < thres):
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in cols]))


if __name__ == '__main__':
    if len(sys.argv) > 3:
        calc_prob(sys.argv[1], sys.argv[2], float(sys.argv[3]))
    else:
        calc_prob(sys.argv[1], sys.argv[2], 1)
