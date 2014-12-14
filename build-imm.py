#! /usr/bin/env python

import sys
import math
import scipy.stats
import numpy as np
import math
import compareResults
import time

debug = False

start = ['ATG']
stop = ['TAA', 'TAG', 'TGA']
nts = ['A', 'T', 'C', 'G']

revCompStart = ['CAT']
revCompStop = ['TTA', 'CTA', 'TCA']

scoreThreshold = -1.30

maxLength = 9
numOrfs = 7

orfMinLength = 90

# If a kmer has at least this many occurences, don't use a smaller kmer
countThreshold = 400

def init():
    counts = []
    for o in xrange(numOrfs):
        countsCurr = []
        for l in xrange(maxLength+1):
            countsCurr.append(dict())
        counts.append(countsCurr)

    # Save imm scores once we've calculated them
    immScores = []
    for i in xrange(numOrfs):
        immScores.append([])
        for j in xrange(maxLength+1):
            immScores[i].append(dict())

    return counts, immScores

def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def findLongORFs(genome):
    minLen = 500

    counts = [0,0]

    longORFs = []
    startPos = None
    for i in xrange(3):
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue

            if startPos == None:# and codon in start:
                startPos = nt
            elif not startPos == None and codon in stop:
                if nt+3-startPos >= minLen:
                    longORFs += [(startPos, nt+3, 0)]
                    #orfs += [(genome[startPos : nt+3], 0)]
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

            if startPos == None:# and codon in revCompStart:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                if startPos-nt+3 >= minLen:
                    longORFs += [(nt-3, startPos, 3)]
                    #orfs += [(genome[nt-3 : startPos], 3)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    longORFs.sort()
    orfs = []
    if longORFs[0][1] <= longORFs[1][0]:
        orfs += [(genome[longORFs[0][0] : longORFs[0][1]], longORFs[0][2])]
        counts[longORFs[0][2]/3] += 1
    for i in xrange(1, len(longORFs)-1):
        if longORFs[i][0] >= longORFs[i-1][1] and longORFs[i][1] <= longORFs[i+1][0]:
            orfs += [(genome[longORFs[i][0] : longORFs[i][1]], longORFs[i][2])]
            counts[longORFs[i][2]/3] += 1
    if longORFs[-2][1] <= longORFs[-1][0]:
        orfs += [(genome[longORFs[-1][0] : longORFs[-1][1]], longORFs[-1][2])]
        counts[longORFs[-1][2]/3] += 1

    # find noncoding segments
    startPos = None
    numEnds = 0
    count = 0
    noncoding = []
    for nt in xrange(len(genome)):
        codon = genome[nt:nt+3]
        if len(codon) < 3:
            continue

        if startPos == None and codon in stop:
            startPos = nt
            numEnds = 0
        elif not startPos == None and codon in start:
            if numEnds == 0:
                numEnds += 1
            else:
                if nt-startPos >= minLen:
                    noncoding += [genome[startPos:nt]]
                    count += 1
                startPos = None
                numEnds = 0
        elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
            startPos = None
            numEnds = 0

    return orfs, noncoding

def updateORFs(orfs):
    minLen = 500
    orfs = sorted(orfs)
    noncoding = []
    for i in xrange(len(orfs)-1):
        if orfs[i+1][0] - orfs[i][1] >= minLen:
            noncoding += [genome[orfs[i][1]:orfs[i+1][0]]]

    coding = []
    for c in orfs:
        coding += [(genome[c[0]:c[1]], c[2])]

    return coding, noncoding

'''
class orfSet:
    def __init__(self, starts, ends, rf):
        self.starts = starts
        self.ends = ends
        self.rf = rf

    def findBestORF(self):
        
        s = self.starts[0]
        e = self.ends[-1]
        rfScores = score((s,e))
        length = e - s

        if self.rf < 3 and max(rfScores) == rfScores[0]:
            normScore = rfScores[0]/length 
            if normScore > scoreThreshold:
                return (s, e, self.rf, normScore)
        elif self.rf >= 3 and max(rfScores) == rfScores[3]:
            normScore = rfScores[3]/length
            if normScore > scoreThreshold:
                return (s, e, self.rf, normScore)
        return None
        

        best = None
        bestScore = scoreThreshold

        for s in self.starts:
            for e in self.ends:
                rfScores = score((s,e))
                length = e - s

                if self.rf < 3 and max(rfScores) == rfScores[0]:
                    normScore = rfScores[0] / length
                    if normScore > bestScore:
                        best = (s,e)
                        bestScore = normScore
                elif self.rf >= 3 and max(rfScores) == rfScores[3]:
                    normScore = rfScores[3] / length
                    if normScore > bestScore:
                        best = (s,e)
                        bestScore = normScore
        
        #if (len(self.starts) > 1 or len(self.ends) > 1) and not best == None:
        #    print 'Scoring: ' + str(self.starts) + ' & ' + str(self.ends)
        #    print best
        #    print bestScore
        #    print ''

        if best == None:
            return None
        else:
            return (best[0], best[1], self.rf, bestScore)
'''


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

def buildIMM(counts, codingSeqs, noncodingSeqs):
    # frame = 0 (forward) or 3 (reverse complement)
    for (seq,frame) in codingSeqs:
        rf = frame
        for i in xrange(len(seq)):
            for l in xrange(min(maxLength+1, i+1)):
                kmer = seq[i-l:i+1]
                if kmer in counts[rf][l]:
                    counts[rf][l][kmer] += 1
                else:
                    counts[rf][l][kmer] = 1
            rf = (rf+1) % 3 + frame

    rf = 6
    for seq in noncodingSeqs:
        for i in xrange(len(seq)):
            for l in xrange(min(maxLength+1, i+1)):
                kmer = seq[i-l:i+1]
                if kmer in counts[rf][l]:
                    counts[rf][l][kmer] += 1
                else:
                    counts[rf][l][kmer] = 1

    return counts


def buildMarkovChain(counts, codingSeqs, noncodingSeqs, length):
    # frame = 0 (forward) or 3 (reverse complement)
    for (seq,frame) in codingSeqs:
        rf = (length % 3) + frame
        for i in xrange(length, len(seq)):
            kmer = seq[i-length:i+1]
            if kmer in counts[rf][length]:
                counts[rf][length][kmer] += 1
            else:
                counts[rf][length][kmer] = 1
            rf = (rf+1) % 3 + frame

    rf = 6
    for seq in noncodingSeqs:
        for i in xrange(length, len(seq)):
            kmer = seq[i-length:i+1]
            if kmer in counts[rf][length]:
                counts[rf][length][kmer] += 1
            else:
                counts[rf][length][kmer] = 1

    return counts
    

def glimmer(genome, counts, immScores, mcLength=None):
    '''
        mcLength=None for IMM, mcLength=x for Markov Chain
    '''

    # Find orfs for each reading frame
    orfs = []
    for i in xrange(6):
        orfs.append([])

    for i in xrange(3):
        startPos = None
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue
                
            if codon in start and startPos == None:
                startPos = nt
            elif not startPos == None and codon in stop:
                if nt+3-startPos >= orfMinLength:
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

            if codon in revCompStart and startPos == None:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                if startPos-nt+3 >= orfMinLength:
                    orfs[i+3] += [(nt-3, startPos)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    '''
    for i in xrange(3):
        starts = []
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue
                
            if codon in start:
                starts.append(nt)
            elif len(starts) > 0 and codon in stop:
                #orfs[i] += [(starts[0], nt+3)]
                orfs[i].append(orfSet(starts, [nt+3], i))
                starts = []
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                starts = []

    for i in xrange(3):
        starts = []
        for pos in xrange(i, len(genome), 3):
            nt = len(genome)-pos
            codon = genome[nt-3:nt]
            if len(codon) < 3:
                continue

            if codon in revCompStart:
                starts.append(nt)
            elif len(starts) > 0 and codon in revCompStop:
                #orfs[i+3] += [(nt-3, starts[0])]
                orfs[i+3].append(orfSet([nt-3], starts, i+3))
                starts = []
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                starts = []
    '''

    # filter out bad-scoring orfs
    goodORFs = []

    totals = [0]*6
    approved = [0]*6

    for i in xrange(len(orfs)):
        if debug:
            print 'Processing rf %d (%d)' % (i, len(orfs[i]))
        
        for orf in orfs[i]:
            if mcLength == None:
                rfScores = scoreHMM(orf, counts, immScores)
            else:
                rfScores = scoreMarkovChain(orf, mcLength, counts)

            # Find largest and second largest scores
            maxScore = max(rfScores)

            length = orf[1]-orf[0]
            if i < 3 and maxScore == rfScores[0] and rfScores[0]/length > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[0]/length))
                approved[i] += 1
            elif i >= 3 and maxScore == rfScores[3] and rfScores[3]/length > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[3]/length))
                approved[i] += 1
            
            '''
            bounds = orf.findBestORF()
            if not bounds == None:
                goodORFs.append(bounds)
            '''

            totals[i] += 1

    # resolve overlapping orfs
    minOverlap = 20
    goodORFs = sorted(goodORFs)
    approvedORFs = []
    suspectORFs = []

    i = 0
    while i < len(goodORFs)-1:
        j = i+1
        while j < len(goodORFs) and (goodORFs[i][1] - goodORFs[j][0] >= minOverlap):
            if mcLength == None:
                overlapScores = scoreHMM((goodORFs[j][0], goodORFs[i][1]), counts, immScores)
            else:
                overlapScores = scoreMarkovChain((goodORFs[j][0], goodORFs[i][1]), mcLength, counts)

            maxScore = max(overlapScores)

            if (goodORFs[i][1]-goodORFs[i][0]) > (goodORFs[j][1]-goodORFs[j][0]) and maxScore == overlapScores[goodORFs[i][2]]:
                del goodORFs[j]
                j -= 1
            elif (goodORFs[i][1]-goodORFs[i][0]) < (goodORFs[j][1]-goodORFs[j][0]) and maxScore == overlapScores[goodORFs[j][2]]:
                # reject shorter ORF
                del goodORFs[i]
                i -= 1
                break
            j += 1
        i += 1

    suspected = [0]*len(goodORFs)
    for i in xrange(len(goodORFs)):
        j = 1
        while i+j < len(goodORFs) and goodORFs[i][1] - goodORFs[i+j][0] >= minOverlap:
            suspected[i] = 1
            suspected[i+j] = 1
            j += 1

        if suspected[i] == 0:
            approvedORFs.append((goodORFs[i][0], goodORFs[i][1], goodORFs[i][2]))
        else:
            suspectORFs.append((goodORFs[i][0], goodORFs[i][1], goodORFs[i][2]))

    return approvedORFs, suspectORFs

def writeORFs(approvedORFs, suspectORFs):
    print '  %d approved, %d suspected ORFs found' % (len(approvedORFs), len(suspectORFs))

    with open('predicted.txt', 'w') as f:
        for orf in approvedORFs:
            f.write(str(orf[0]) + '\t' + str(orf[1]) + '\n')
        for orf in suspectORFs:
            f.write('*' + str(orf[0]) + '\t' + str(orf[1]) + '\n')

def scoreHMM(orf, counts, immScores):
    # Score the given region in first 3 reading frames
    rfScores = []
    for rf in xrange(3):
        score = 0
        for x in xrange(orf[0], orf[1]):
            score += math.log(imm(rf, genome[max(orf[0], x-maxLength) : x+1], counts, immScores))
            rf = (rf+1) % 3
        rfScores.append(score)

    # Score the given region in reverse 3 reading frames
    for rf in xrange(3):
        score = 0
        for x in xrange(orf[0], orf[1]):
            score += math.log(imm(rf+3, genome[max(orf[0], x-maxLength) : x+1], counts, immScores))
            rf = (rf+1) % 3
        rfScores.append(score)

    
    # Score in last reading frame
    score = 0
    for x in xrange(orf[0], orf[1]):
        score += math.log(imm(6, genome[max(orf[0], x-maxLength) : x+1], counts, immScores))
    rfScores.append(score)

    return rfScores

def scoreMarkovChain(orf, length, counts):
    # Score the given region in first 3 reading frames
    rfScores = []
    for rf in xrange(3):
        score = 0
        rf = (rf+length) % 3
        for x in xrange(orf[0]+length, orf[1]):
            score += math.log(prob(rf, genome[x-length : x+1], counts))
            rf = (rf+1) % 3
        rfScores.append(score)

    # Score the given region in reverse 3 reading frames
    for rf in xrange(3):
        score = 0
        rf = (rf+length) % 3
        for x in xrange(orf[0], orf[1]):
            score += math.log(prob(rf+3, genome[x-length : x+1], counts))
            rf = (rf+1) % 3
        rfScores.append(score)

    
    # Score in last reading frame
    score = 0
    for x in xrange(orf[0]+length, orf[1]):
        score += math.log(prob(6, genome[x-length : x+1], counts))
    rfScores.append(score)

    return rfScores

def imm(rf, seq, counts, immScores):
    if seq in immScores[rf][len(seq)-1]:
        return immScores[rf][len(seq)-1][seq]

    prefix = seq[:-1]

    if len(seq) == 1:
        probs = [0]*4
        for i in xrange(4):
            immScores[rf][0][prefix+nts[i]] = prob(rf, prefix+nts[i], counts)
        return immScores[rf][len(seq)-1][seq]

    length = len(seq)-1
    #if seq in counts[rf][length] and counts[rf][length][seq] > countThreshold:
    if prefix in counts[rf][length-1] and counts[rf][length-1][prefix] > countThreshold:
        wgt = 1
    else:
        '''
        currCounts = [1]*4
        nextCounts = [1]*4
        for i in xrange(4):
            tempSeq = prefix + nts[i]
            if tempSeq in counts[rf][length]:
                currCounts[i] += counts[rf][length][tempSeq]

            tempSeq = prefix[1:] + nts[i]
            if tempSeq in counts[rf][length-1]:
                nextCounts[i] += counts[rf][length-1][tempSeq]
        '''

        currCounts = [1]*4
        currTotal = 4
        for i in xrange(4):
            tempSeq = prefix + nts[i]
            if tempSeq in counts[rf][length]:
                currCounts[i] += counts[rf][length][tempSeq]
                currTotal += counts[rf][length][tempSeq]

        #currProbs = [0]*4
        nextProbs = [0]*4
        for i in xrange(4):
            #currProbs[i] = float(currCounts[i]) / float(currTotal)
            nextProbs[i] = imm(rf, prefix[1:]+nts[i], counts, immScores)
        
        #d = 1 - scipy.stats.chi2_contingency(np.array([currCounts, nextCounts]))[1]
        d = scipy.stats.chi2_contingency(np.array([currCounts, nextProbs]))[1]

        if d < 0.5:
            wgt = 0
        else:
            if prefix in counts[rf][length-1]:
                wgt = d * counts[rf][length-1][prefix] / 400.0
            else:
                wgt = 0

    #score = wgt * prob(rf, seq, counts) + (1-wgt) * imm(rf, seq[1:], counts, immScores)
    probs = [0]*4
    for i in xrange(4):
        probs[i] = wgt * prob(rf, prefix+nts[i], counts) + (1-wgt) * imm(rf, prefix[1:]+nts[i], counts, immScores)
    totalProb = sum(probs)

    for i in xrange(4):
        immScores[rf][len(seq)-1][prefix+nts[i]] = float(probs[i]) / totalProb

    return immScores[rf][len(seq)-1][seq]

def prob(rf, seq, counts):
    length = len(seq)-1

    # Count all possible next nucleotides, with +1 smoothing
    totalNum = 4
    for n in nts:
        if seq[:-1]+n in counts[rf][length]:
            totalNum += counts[rf][length][seq[:-1]+n]

    if not seq in counts[rf][length]:
        return 1.0 / float(totalNum)
    else:
        return float(1+counts[rf][length][seq]) / float(totalNum)

def runIMM():
    print 'IMM max length %d:' % maxLength
    #maxLength = maxImmLength
    counts, immScores = init()

    startTime = time.time()
    coding, noncoding = findLongORFs(genome)
    counts = buildIMM(counts, coding, noncoding)
    buildTime = time.time()
    approved, suspect = glimmer(genome, counts, immScores)
    endTime = time.time()
    print '  %0.2fs total (%0.2fs to build, %0.2fs to run)' % (endTime-startTime, buildTime-startTime, endTime-buildTime)
    writeORFs(approved, suspect)

    compareResults.compare('trueORFs.txt', 'predicted.txt')
    print ''

def runMC(length):
    print 'Markov Chain length %d:' % length
    counts, immScores = init()

    startTime = time.time()
    coding, noncoding = findLongORFs(genome)
    counts = buildMarkovChain(counts, coding, noncoding, length)
    buildTime = time.time()
    approved, suspect = glimmer(genome, counts, immScores, length)
    endTime = time.time()
    print '  %0.2fs total (%0.2fs to build, %0.2fs to run)' % (endTime-startTime, buildTime-startTime, endTime-buildTime)
    writeORFs(approved, suspect)

    compareResults.compare('trueORFs.txt', 'predicted.txt')
    print ''

def runIterativeIMM():
    ''' Repeatedly train on the set of predicted ORFs, and use the IMM predict a new set of ORFs.
        Continue until the set of predicted ORFs does not change.
    '''

    print 'Iterative IMM max length %d:' % maxLength

    coding, noncoding = findLongORFs(genome)
    oldORFs = []
    i = 0
    for x in xrange(10):
        i += 1
        print '  Iteration %d' % i

        counts, immScores = init()
        startTime = time.time()
        print '    Training on %d coding, %d noncoding...' % (len(coding), len(noncoding))
        counts = buildIMM(counts, coding, noncoding)
        buildTime = time.time()
        print '      Done, %0.2fs' % (buildTime-startTime)
        print '    Recalculating ORFs...'
        approved, suspect = glimmer(genome, counts, immScores)
        endTime = time.time()
        print '      Done, %0.2fs' % (endTime-buildTime)
        print '    Found %d approved, %d suspected' % (len(approved), len(suspect))

        writeORFs(approved, suspect)
        compareResults.compare('trueORFs.txt', 'predicted.txt')
        print ''

        newORFs = approved+suspect
        if noChange(oldORFs, newORFs):
            print '  Converged'
            break
        else:
            coding, noncoding = updateORFs(approved)
            oldORFs = newORFs
    print 'Completed after %d iterations' % i

def runIterativeMC(length):
    ''' Repeatedly train on the set of predicted ORFs, and use a fixed-length MM predict a new set of ORFs.
        Continue until the set of predicted ORFs does not change.
    '''

    print 'Iterative MC length %d:' % length

    coding, noncoding = findLongORFs(genome)
    oldORFs = []
    i = 0
    for x in xrange(10):
        i += 1
        print '  Iteration %d' % i

        counts, immScores = init()
        startTime = time.time()
        print '    Training on %d coding, %d noncoding...' % (len(coding), len(noncoding))
        counts = buildMarkovChain(counts, coding, noncoding, length)
        buildTime = time.time()
        print '      Done, %0.2fs' % (buildTime-startTime)
        print '    Calculating ORFs...'
        approved, suspect = glimmer(genome, counts, immScores, length)
        endTime = time.time()
        print '      Done, %0.2fs' % (endTime-buildTime)
        print '    Found %d approved, %d suspected' % (len(approved), len(suspect))

        writeORFs(approved, suspect)
        compareResults.compare('trueORFs.txt', 'predicted.txt')
        print ''

        newORFs = approved+suspect
        if noChange(oldORFs, newORFs):
            print '  Converged'
            break
        else:
            coding, noncoding = updateORFs(approved)
            oldORFs = newORFs
    print 'Completed after %d iterations' % i

def noChange(oldORFs, newORFs):
    if not len(oldORFs) == len(newORFs):
        return False

    oldORFs.sort()
    newORFs.sort()
    return (oldORFs == newORFs)

genome = readFASTA(sys.argv[1])


#for i in xrange(maxLength+1):
#    runMC(i)
#for i in xrange(maxLength+1):
#    maxLength = i
#    runIMM()

'''
for i in xrange(8):
    maxLength = 8
    scoreThreshold = -1.3 + 0.01*i
    print 'scoreThreshold = %0.2f' % (scoreThreshold)
    runIMM() 
'''

maxLength = 8
runIMM()
#runMC(6)

#runIterativeIMM()
#runIterativeMC(6)





    