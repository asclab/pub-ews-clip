#!/usr/bin/env python
#!/usr/bin/env python
import sys

genes = {}

name = ''
alpha = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for line in sys.stdin:
    spl = line.strip().split('=')
    if spl[0] == 'SEQUENCE_ID':
        name = spl[1].split(';')[0]
        if not name in genes:
            genes[name] = 0
        else:
            genes[name] += 1

    if spl[0] == 'PRIMER_LEFT_0_SEQUENCE':
        sys.stdout.write('%s_%s_fwd\t%s\n' % (name, alpha[genes[name]], spl[1]))
    if spl[0] == 'PRIMER_RIGHT_0_SEQUENCE':
        sys.stdout.write('%s_%s_rev\t%s\n' % (name, alpha[genes[name]], spl[1]))
