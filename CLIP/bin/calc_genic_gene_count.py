#!/usr/bin/env python
#!/usr/bin/env python
'''
Calculate the number of reads for each gene for EWS1 and EWS2 and HeLa RNA by genic region
'''

import sys
import gzip


def write(genenames, genecount_ews1, genecount_ews2, genecount_hela):
    for k in sorted(genenames):
        if genecount_ews1[k] + genecount_ews2[k] == 0:
            continue
        cols = [k[0], k[1], genenames[k], genecount_ews1[k], genecount_ews2[k], genecount_hela[k]]
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in cols]))




def openfile(fname):
    if fname[-3:] == '.gz':
        return gzip.open(fname)
    return open(fname)


def calc_genecount(fname, thres=0.01):
    ref = None
    genenames = {}
    genecount_hela = {}
    genecount_ews1 = {}
    genecount_ews2 = {}

    ews1_col = 4
    ews2_col = 5
    hela_col = 6
    gene_col = 8
    genic_region_col = 9
    gene_name_col = 10

    header = True
    sys.stderr.write("Calculating gene totals\n")
    i = 0

    with openfile(fname) as f:
        for line in f:
            i+=1
            if not line.strip() or line[0] == '#':
                continue

            if header:
                header = False
                continue

            cols = line.strip('\n').split('\t')

            if ref and cols[0] != ref:
                write(genenames, genecount_ews1, genecount_ews2, genecount_hela)
                genenames = {}
                genecount_hela = {}
                genecount_ews1 = {}
                genecount_ews2 = {}

            if cols[0] == 'chrM':
                continue

            ref = cols[0]

            if i > 100000:
                sys.stderr.write('%s:%s\n' % (cols[0], cols[1]))
                i = 0

            if not cols[gene_col]:
                continue
        
            genes = cols[gene_col].split(',')
            names = cols[gene_name_col].split(',')
            region = cols[genic_region_col]

            ews1_count = int(cols[ews1_col]) if cols[ews1_col] else 0
            ews2_count = int(cols[ews2_col]) if cols[ews2_col] else 0
            hela_count = int(cols[hela_col]) if cols[hela_col] else 0
            
            for gene, name in zip(genes, names):
                if not (gene,region) in genenames:
                    genenames[(gene,region)] = name
                    genecount_ews1[(gene,region)] = 0
                    genecount_ews2[(gene,region)] = 0
                    genecount_hela[(gene,region)] = 0

                genecount_ews1[(gene,region)] += ews1_count
                genecount_ews2[(gene,region)] += ews2_count
                genecount_hela[(gene,region)] += hela_count

    if genenames:
        write(genenames, genecount_ews1, genecount_ews2, genecount_hela)




if __name__ == '__main__':
    calc_genecount(sys.argv[1])

