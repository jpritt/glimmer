#! /usr/bin/env python

import sys

def readORFs(filename):
    ''' Given a list of coding sequences from GenBank, writes a list of genes with exons clearly marked to an easily readable file.
        This method is intended for use with build-imm.py and is primarily meant for prokaryotic genomes, as it does not retain any
         intron information.
        The input file should be an annotated file containing coding sequences e.g. as available from GenBank.
    '''
    orfs = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                # All ORFs begin with '>'
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

                # Read the start and end positions of the ORF from the location string
                if locTagEnd > locTagStart:
                    value = line[locTagStart:locTagEnd]

                    startLocA = value.find('..')
                    startLocB = value.find(',')
                    if startLocB < 0 or (startLocB > startLocA):
                        start = value[:startLocA]
                    elif startLocA < 0 or (startLocB < startLocA):
                        start = value[:startLocB]
                    else:
                        print 'Error! Line:'
                        print value
                        exit()

                    while not start[-1].isdigit():
                        start = start[:-1]

                    endLocA = value.rfind('..')
                    endLocB = value.rfind(',')
                    end = value[max(endLocA, endLocB)+2:]
                    while not end[0].isdigit():
                        end = end[1:]
                    #if join:
                    #    end = value[value.rfind('..')+2:]
                    #else:
                    #    end = value[value.index('..')+2:]


                    # Write this ORF as a single line
                    if complement:
                        orfs.append((int(start)-1,end))
                    else:
                        orfs.append((int(start)+1,end))
    return orfs

orfs = readORFs(sys.argv[1])
with open(sys.argv[2], 'w') as f:
    for (a,b) in orfs:
        f.write(str(int(a)) + '\t' + str(b) + '\n')

