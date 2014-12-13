#! /usr/bin/env python

import sys

def readORFs(filename):
    orfs = []
    with open(filename, 'r') as f:
        for line in f:
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

                if locTagEnd > locTagStart:
                    value = line[locTagStart:locTagEnd]
                    start = value[:value.index('..')]
                    end = value[value.index('..')+2:]

                    if complement:
                        orfs.append((int(start)-1,end))
                    else:
                        orfs.append((int(start)+1,end))
    return orfs

orfs = readORFs(sys.argv[1])
with open('trueORFs.txt', 'w') as f:
    for (a,b) in orfs:
        f.write(str(int(a)) + '\t' + str(b) + '\n')

