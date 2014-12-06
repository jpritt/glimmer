#! /usr/bin/env python

import sys
with open(sys.argv[1], 'r') as f:
    with open('sequenceSmaller.fasta', 'w') as f2:
        i = 0
        for line in f:
            f2.write(line)

            i += 1
            if i > 5000:
                break

