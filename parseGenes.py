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

                join = False
                if line[locTagStart:locTagStart+5] == 'join(':
                    join = True
                    locTagStart += 5
                    locTagEnd -= 1


                while not line[locTagStart].isdigit():
                    locTagStart += 1
                while not line[locTagEnd-1].isdigit():
                    locTagEnd -= 1

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

                    if complement:
                        orfs.append((int(start)-1,end))
                    else:
                        orfs.append((int(start)+1,end))
    return orfs

orfs = readORFs(sys.argv[1])
with open('trueORFs.txt', 'w') as f:
    for (a,b) in orfs:
        f.write(str(int(a)) + '\t' + str(b) + '\n')

