#! /usr/bin/env python

import sys

with open(sys.argv[1], 'r') as f:
    domains = []
    for line in f:
        row = line.rstrip().split('\t')

        if len(row) < 5:
            continue

        domainType = row[2]
        start = int(row[3])
        end = int(row[4])
        domains.append((start, end, domainType))

total = 0
correct = 0
with open(sys.argv[2], 'r') as f:
    for line in f:
        row = line.rstrip().split('\t')

        if not len(row) == 2:
            continue

        start = int(row[0])
        end = int(row[1])

        for (a,b,c) in domains:
            if start == a and end == b:
                correct += 1
                print 'Found %d-%d (%s)' % (a,b,c)
                break
        total += 1

print '%d / %d correct' % (correct, total)
