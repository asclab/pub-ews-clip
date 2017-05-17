#!/usr/bin/env python

import sys

lo = int(sys.argv[2])
hi = int(sys.argv[3])

class_count = 0
outofbounds = 0

bait = sys.argv[1].split(',')

passed = {}

for b in bait:
    passed[b] = set()

failed = set()

with open(sys.argv[4]) as f:
    for line in f:
        cols = line.strip().split('\t')
        if cols[1] in bait:
            class_count += 1
            if lo <= int(cols[2]) <= hi:
                passed[cols[1]].add(cols[0])
            else:
                outofbounds += 1
        else:
            failed.add(cols[0])

sys.stderr.write('%s: %s\n' % (sys.argv[1], class_count))
for b in bait:
    sys.stderr.write('Passed [%s]: %s\n' % (b, len(passed[b])))
sys.stderr.write('OOB: %s\n' % (outofbounds))
sys.stderr.write('Failed: %s\n' % len(failed))

######
# check for membership in all classes - not needed...
######
# count = 0
# with open(sys.argv[5]) as f:
#     for line in f:
#         if line[0] == '#':
#             continue
#         cols = line.strip().split('\t')
#         found = True
#         for b in passed:
#             if not cols[3] in passed[b]:
#                 found = False

#         if found and cols[3] not in failed:
#             sys.stdout.write(line)
#             count += 1
#             #sys.stdout.write('%s\t%s\t%s\n' % (cols[0], cols[1], cols[2]))

count = 0
with open(sys.argv[5]) as f:
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')

        found = False
        for b in bait:
            if cols[3] in passed[b]:
                found = True

        if found:
            sys.stdout.write(line)
            count += 1


sys.stderr.write('Written: %s\n' % count)

