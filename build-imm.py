#! /usr/bin/env python

import sys
import math
import scipy.stats
import numpy as np

start = ['ATG']
stop = ['TAA', 'TAG', 'TGA']
nts = ['A', 'T', 'C', 'G']

scoreThreshold = -1000

maxLength = 8
numOrfs = 6
counts = []
for o in xrange(numOrfs):
    countsCurr = []
    for l in xrange(maxLength+1):
        countsCurr.append(dict())
    counts.append(countsCurr)

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

    for i in xrange(3):
        startPos = None
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3] 
            if codon in start:
                startPos = nt
            elif codon in stop and not startPos == None:
                if nt-startPos >= minLen:
                    orfs += [(genome[startPos : nt], 0)]
                startPos = None

    rev = reverseComplement(genome)
    for i in xrange(3):
        startPos = None
        for nt in xrange(i, len(rev), 3):
            codon = rev[nt:nt+3] 
            if codon in start:
                startPos = nt
            elif codon in stop and not startPos == None:
                if nt-startPos >= minLen:
                    orfs += [(genome[startPos : nt], 0)]
                startPos = None

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
            print 'Unrecognized nucleotide: ' + c
            exit()
    return rev

def buildIMM(seqs):
    for (seq,rf) in seqs:
        for i in xrange(len(seq)):
            for l in xrange(min(maxLength+1, i+1)):
                kmer = seq[i-l:i+1]
                if kmer in counts[rf][l]:
                    counts[rf][l][kmer] += 1
                else:
                    counts[rf][l][kmer] = 1
            rf = (rf+1) % 3

    print len(counts)
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

    for i in xrange(3):
        starts = []
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3] 
            if codon in start:
                starts += [nt]
            elif codon in stop and len(starts) > 0:
                orfs[i] += [(starts[0], nt)]
                starts = []

    rev = reverseComplement(genome)
    for i in xrange(3):
        starts = []
        for nt in xrange(i, len(rev), 3):
            codon = rev[nt:nt+3] 
            if codon in start:
                starts += [nt]
            elif codon in stop and len(starts) > 0:
                orfs[i+3] += [(starts[0], nt)]
                starts = []

    # filter out bad-scoring orfs
    goodORFs = []
    for i in xrange(len(orfs)):
        for orf in orfs[i]:
            rfScores = score(orf)
            if rfScores[0] > rfScores[1] and rfScores[0] > rfScores[2] and rfScores[0] > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[0]))


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
    # Score the given region in all 3 reading frames
    rfScores = []
    for rf in xrange(3):
        score = 0
        for x in xrange(orf[0], orf[1]):
            score += math.log(imm(rf, genome[max(orf[0], x-maxLength+1) : x+1]))
            rf = (rf+1) % 3
        rfScores.append(score)
    return rfScores

def imm(rf, seq):
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

    return wgt * prob(rf, seq) + (1-wgt) * imm(rf, seq[1:])

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




    