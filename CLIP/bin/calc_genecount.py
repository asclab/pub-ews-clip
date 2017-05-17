#!/usr/bin/env python
#!/usr/bin/env python
'''
Calculate the number of reads for each gene for EWS1 and EWS2

Ignore anti-sense bins!

'''

import sys
import gzip


def openfile(fname):
    if fname[-3:] == '.gz':
        return gzip.open(fname)
    return open(fname)


def calc_genecount(fname, thres=0.01):
    genenames = {}
    genecount_ews1 = {}
    genecount_ews2 = {}

    ews1_col = 4
    ews2_col = 5
    gene_col = 8
    region_col = 9
    gene_name_col = 10

    header = True
    sys.stderr.write("Calculating gene totals\n")
    i = 0

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
            regions = cols[region_col].split(',')

            is_anti = True
            for reg in regions:
                if '_ANTI' not in reg:
                    is_anti = False

            if is_anti:
                continue

            ews1_count = int(cols[ews1_col]) if cols[ews1_col] else 0
            ews2_count = int(cols[ews2_col]) if cols[ews2_col] else 0

            for gene, name in zip(genes, names):
                if not gene in genenames:
                    genenames[gene] = name
                    genecount_ews1[gene] = 0
                    genecount_ews2[gene] = 0

                genecount_ews1[gene] += ews1_count
                genecount_ews2[gene] += ews2_count

    for gene in genenames:
        if genecount_ews1[gene] + genecount_ews2[gene] == 0:
            continue
        cols = [gene, genenames[gene], genecount_ews1[gene], genecount_ews2[gene]]
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in cols]))


if __name__ == '__main__':
    calc_genecount(sys.argv[1])
