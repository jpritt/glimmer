#! /usr/bin/env python

import sys
import math
import scipy.stats
import numpy as np

start = ['ATG']
stop = ['TAA', 'TAG', 'TGA']
nts = ['A', 'T', 'C', 'G']

revCompStart = ['CAT']
revCompStop = ['TTA', 'CTA', 'TCA']

scoreThreshold = -10000

maxLength = 8
numOrfs = 6
counts = []
for o in xrange(numOrfs):
    countsCurr = []
    for l in xrange(maxLength+1):
        countsCurr.append(dict())
    counts.append(countsCurr)

orfMinLength = 50

# Save imm scores once we've calculated them
immScores = []
for i in xrange(6):
    immScores.append([])
    for j in xrange(maxLength):
        immScores[i].append(dict())

# If a kmer has at least this many occurences, don't use a smaller kmer
countThreshold = 400

def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def findLongORFs(genome):
    minLen = 500
    orfs = []

    counts = [0,0]

    startPos = None
    for i in xrange(3):
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue

            if startPos == None and codon in start:
                startPos = nt
            elif not startPos == None and codon in stop:
                if nt+3-startPos >= minLen:
                    orfs += [(genome[startPos : nt+3], 0)]
                    counts[0] += 1
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    for i in xrange(3):
        startPos = None
        for pos in xrange(i, len(genome), 3):
            nt = len(genome)-pos
            codon = genome[nt-3:nt]
            if len(codon) < 3:
                continue

            if startPos == None and codon in revCompStart:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                if startPos-nt+3 >= minLen:
                    orfs += [(genome[nt-3 : startPos], 3)]
                    counts[1] += 1
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    print counts
    return orfs

def reverseComplement(seq):
    rev = ''
    for c in seq:
        if c == 'G':
            rev = 'C' + rev
        elif c == 'C':
            rev = 'G' + rev
        elif c == 'A':
            rev = 'T' + rev
        elif c == 'T':
            rev = 'A' + rev
        else:
            rev = 'N' + rev
            #print 'Unrecognized nucleotide: ' + c
            #exit()
    return rev

def buildIMM(seqs):
    # frame = 0 (forward) or 3 (reverse complement)
    for (seq,frame) in seqs:
        rf = frame
        for i in xrange(len(seq)):
            for l in xrange(min(maxLength+1, i+1)):
                kmer = seq[i-l:i+1]
                if kmer in counts[rf][l]:
                    counts[rf][l][kmer] += 1
                else:
                    counts[rf][l][kmer] = 1
            rf = (rf+1) % 3 + frame

    print counts[0][0]
    print counts[1][0]
    print counts[2][0]
    print counts[3][0]
    print counts[4][0]
    print counts[5][0]

def glimmer(genome):
    # Find orfs for each reading frame
    orfs = []
    for i in xrange(6):
        orfs.append([])

    startPos = None
    for i in xrange(3):
        startPos = None
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue
                
            if startPos == None and codon in start:
                startPos = nt
            elif not startPos == None and codon in stop:
                orfs[i] += [(startPos, nt+3)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    for i in xrange(3):
        startPos = None
        for pos in xrange(i, len(genome), 3):
            nt = len(genome)-pos
            codon = genome[nt-3:nt]
            if len(codon) < 3:
                continue

            if startPos == None and codon in revCompStart:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                orfs[i+3] += [(nt-3, startPos)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    #with open('predicted.txt', 'w') as f:
    #    for i in xrange(6):
    #        for orf in orfs[i]:
    #            f.write(str(orf[0]) + '\t' + str(orf[1]) + '\n')
    #exit()

    # filter out bad-scoring orfs
    goodORFs = []

    totals = [0]*6
    approved = [0]*6

    for i in xrange(len(orfs)):
        print 'Processing rf %d (%d)' % (i, len(orfs[i]))
        
        for orf in orfs[i]:
            rfScores = score(orf)
            if i < 3 and max(rfScores) == rfScores[0] and rfScores[0] > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[0]))
                approved[i] += 1
            elif i >= 3 and max(rfScores) == rfScores[3] and rfScores[3] > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[0]))
                approved[i] += 1

            totals[i] += 1

    print 'Approved: ' + str(approved)
    print 'Totals:   ' + str(totals)


    print '%d potential regions found' % len(goodORFs)

    # resolve overlapping orfs
    goodORFs = sorted(goodORFs)
    approvedORFs = []
    suspectORFs = []
    i = 0
    while i < len(goodORFs)-1:
        if goodORFs[i][1] > goodORFs[i+1][0]:
            overlapScores = score((goodORFs[i+1][0], goodORFs[i][1])) #score(genome[goodORFs[i+1][0]:genome[goodORFs[i][1]]])
            maxScore = 0
            if overlapScores[1] > overlapScores[0] and overlapScores[1] > overlapScores[2]:
                maxScore = 1
            elif overlapScores[2] > overlapScores[0] and overlapScores[2] > overlapScores[1]:
                maxScore = 2

            if (goodORFs[i][1]-goodORFs[i][0]) > (goodORFs[i+1][1]-goodORFs[i+1][0]) and maxScore == goodORFs[i][2]:
                # reject shorter ORF
                approvedORFs.append((goodORFs[i][0], goodORFs[i][1]))
            elif (goodORFs[i][1]-goodORFs[i][0]) < (goodORFs[i+1][1]-goodORFs[i+1][0]) and maxScore == goodORFs[i+1][2]:
                # reject shorter ORF
                approvedORFs.append((goodORFs[i+1][0], goodORFs[i+1][1]))
            else:
                suspectORFs.append((goodORFs[i][0], goodORFs[i][1]))
                suspectORFs.append((goodORFs[i+1][0], goodORFs[i+1][1]))
            i += 2
        else:
            approvedORFs.append((goodORFs[i][0], goodORFs[i][1]))
            i += 1
    # append the last orf
    approvedORFs.append((goodORFs[-1][0], goodORFs[-1][1]))

    print '%d approved, %d suspected ORFs found' % (len(approvedORFs), len(suspectORFs))

    with open('predicted.txt', 'w') as f:
        for orf in approvedORFs:
            f.write(str(orf[0]) + '\t' + str(orf[1]) + '\n')
        for orf in suspectORFs:
            f.write(str(orf[0]) + '\t' + str(orf[1]) + '\n')

def score(orf):
    # Score the given region in first 3 reading frames
    rfScores = []
    for rf in xrange(3):
        score = 0
        for x in xrange(orf[0], orf[1]):
            score += math.log(imm(rf, genome[max(orf[0], x-maxLength+1) : x+1]))
            rf = (rf+1) % 3
        rfScores.append(score)

    # Score the given region in reverse 3 reading frames
    for rf in xrange(3):
        score = 0
        for x in xrange(orf[0], orf[1]):
            score += math.log(imm(rf+3, genome[max(orf[0], x-maxLength+1) : x+1]))
            rf = (rf+1) % 3
        rfScores.append(score)
    return rfScores

def imm(rf, seq):
    if seq in immScores[rf][len(seq)-1]:
        return immScores[rf][len(seq)-1][seq]

    if len(seq) == 1:
        return prob(rf, seq)

    length = len(seq)-1
    if seq in counts[rf][length] and counts[rf][length][seq] > countThreshold:
        wgt = 1
    else:
        currCounts = [1]*4
        nextCounts = [1]*4
        for i in xrange(4):
            tempSeq = seq[:-1] + nts[i]
            if tempSeq in counts[rf][length]:
                currCounts[i] = counts[rf][length][tempSeq]
            tempSeq = seq[1:-1] + nts[i]
            if tempSeq in counts[rf][length-1]:
                nextCounts[i] = counts[rf][length-1][tempSeq]
        
        c = 1 - scipy.stats.chi2_contingency(np.array([currCounts, nextCounts]))[1]

        if c < 0.5:
            wgt = 0
        else:
            if seq in counts[rf][length]:
                wgt = c * counts[rf][length][seq] / 400.0
            else:
                wgt = 0

    score = wgt * prob(rf, seq) + (1-wgt) * imm(rf, seq[1:])

    immScores[rf][len(seq)-1][seq] = score
    return score

def prob(rf, seq):
    length = len(seq)-1

    if not seq in counts[rf][length]:
        return 0

    totalNum = 0
    for n in nts:
        if seq[:-1]+n in counts[rf][length]:
            totalNum += counts[rf][length][seq[:-1]+n]
    return float(counts[rf][length][seq]) / float(totalNum)


genome = readFASTA(sys.argv[1])
orfs = findLongORFs(genome)
buildIMM(orfs)
glimmer(genome)




    