#! /usr/bin/env python

import sys

'''
    Given a list of coding sequences from GenBank, writes a list of genes with exons clearly marked to an easily readable file
'''
def readORFs(inFile, outFile):
    orfs = []
    with open(outFile, 'w') as f:
        with open(inFile, 'r') as f2:
            for line in f2:
                if line[0] == '>':
                    proteinTagStart = line.index('[protein=')
                    proteinTagEnd = line.index(']', proteinTagStart)
                    if 'hypothetical' in line[proteinTagStart:proteinTagEnd]:
                        continue

                    locTagStart = line.index('[location=') + 10
                    locTagEnd = line.index(']', locTagStart)

                    complement = False
                    if line[locTagStart:locTagStart+11] == 'complement(':
                        complement = True
                        locTagStart += 11
                        locTagEnd -= 1

                    join = False
                    if line[locTagStart:locTagStart+5] == 'join(':
                        join = True
                        locTagStart += 5
                        locTagEnd -= 1


                    while not line[locTagStart].isdigit():
                        locTagStart += 1
                    while not line[locTagEnd-1].isdigit():
                        locTagEnd -= 1

                    segments = line[locTagStart:locTagEnd].split(',')
                    exons = []
                    skip = False
                    for s in segments:
                        bounds = s.split('..')
                        if not len(bounds) == 2:
                            skip = True
                            break

                        for i in xrange(len(bounds)):
                            while not bounds[i][0].isdigit():
                                bounds[i] = bounds[i][1:]
                            while not bounds[i][-1].isdigit():
                                bounds[i] = bounds[i][:-1]

                        exons.append(','.join(bounds))
                    if not skip:
                        f.write('\t'.join(exons) + '\n')

readORFs(sys.argv[1], sys.argv[2])

