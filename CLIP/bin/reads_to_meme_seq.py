#!/usr/bin/env python
#
# Takes a SAM (samtools view) file as stdin, a read list as argv[1], and a fasta file as argv[2]
#
# Return a FASTA file with the read sequence, +/- 50bp, for use in MEME
#

import sys
import subprocess


refname = sys.argv[2]
expand = int(sys.argv[3])

bg = False

if len(sys.argv) > 4:
    if sys.argv[4] == 'bg':
        bg = True

reads = set()
with open(sys.argv[1]) as f:
    for line in f:
        reads.add(line.strip())

def revcomp(seq):
    out = ''
    for base in seq.upper()[::-1]:
        if base == 'A':
            out += 'T'
        elif base == 'T':
            out += 'A'
        elif base == 'C':
            out += 'G'
        elif base == 'G':
            out += 'C'
        elif base == 'N':
            out += 'N'
        else:
            sys.stderr.write('revcomp error: %s\n' % seq)
            sys.exit(1)
    return out


def getseq(chrom, start, end):
    cmd = ['samtools', 'faidx', refname, '%s:%s-%s' % (chrom, start+1, end)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    seq = ""
    while proc.returncode is None:
        out = proc.communicate()[0]
        if out is not None:
            seq += out
        proc.poll()

    ret = ""
    for line in seq.split('\n'):
        if line and line[0] != '>':
            ret += line

    return ret


count = 0
for line in sys.stdin:
    cols = line.strip().split('\t')

    if not bg and cols[0] not in reads:
        continue
    elif bg and cols[0] in reads:
        continue

    count += 1
    flags = int(cols[1])

    chrom = cols[2]
    start = int(cols[3])-1
    end = start
    cigar = []
    clen = ''
    for s in cols[5]:
        if s in '0123456789':
            clen += s
        else:
            cigar.append((s, int(clen)))
            clen = ''

    seq = cols[9]
    if cigar[0][0] == 'S':
        seq = seq[cigar[0][1]:]
    if cigar[-1][0] == 'S':
        seq = seq[:-cigar[-1][1]]



    for op, oplen in cigar:
        if op in 'MND':
            end += oplen

    pre = getseq(chrom, start - expand, start)
    post = getseq(chrom, end, end + expand)

    if flags & 0x10 > 0:
        print ">%s %s:%s-%s" % (cols[0], chrom, end + expand, start - expand+1)
        print revcomp(pre+seq+post)
    else:
        print ">%s %s:%s-%s" % (cols[0], chrom, start - expand+1, end + expand)
        print pre+seq+post

    if count % 10000 == 0:
    	sys.stderr.write('%s:%s %s\n' % (chrom, cols[3], cols[0]))


    #print '\t'.join([str(x) for x in [cols[0], flags,chrom, start+1, end, cigar, cols[9], pre, seq, post, full]])

sys.stderr.write("%s reads\n" % count)