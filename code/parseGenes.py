#! /usr/bin/env python

import sys

def readORFs(inFile, outFile):
    ''' Given a list of coding sequences from GenBank, writes a list of genes with exons clearly marked to an easily readable file.
        This method is intended for eukaryotic genomes for use with HIMM.
        The input file should be an annotated file containing coding sequences e.g. as available from GenBank.
    '''
    orfs = []
    with open(outFile, 'w') as f:
        with open(inFile, 'r') as f2:
            for line in f2:
                # All ORFs begin with '>'
                if line[0] == '>':
                    proteinTagStart = line.index('[protein=')
                    proteinTagEnd = line.index(']', proteinTagStart)

                    # Ignore hypothetical proteins, because these may actually be wrong
                    if 'hypothetical' in line[proteinTagStart:proteinTagEnd]:
                        continue

                    # This tag contains location information for the protein
                    locTagStart = line.index('[location=') + 10
                    locTagEnd = line.index(']', locTagStart)

                    complement = False
                    if line[locTagStart:locTagStart+11] == 'complement(':
                        complement = True
                        locTagStart += 11
                        locTagEnd -= 1

                    # Joined ORFs contain multiple exons
                    join = False
                    if line[locTagStart:locTagStart+5] == 'join(':
                        join = True
                        locTagStart += 5
                        locTagEnd -= 1

                    # Remove any extra symbols
                    while not line[locTagStart].isdigit():
                        locTagStart += 1
                    while not line[locTagEnd-1].isdigit():
                        locTagEnd -= 1

                    # Parse the list of bounds
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

                    # Write this ORF as a single line
                    if not skip:
                        f.write('\t'.join(exons) + '\n')

readORFs(sys.argv[1], sys.argv[2])

