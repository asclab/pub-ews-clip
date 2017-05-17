#!/usr/bin/env python
import sys

name = ''
seq = ''
for line in sys.stdin:
    if line[0] == '>':
        if seq:
            sys.stdout.write('SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=\n' % (name, seq))
        name = line.strip()[1:]
        seq = ''
    else:
        seq += line.strip()

if seq:
    sys.stdout.write('SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=\n' % (name, seq))
