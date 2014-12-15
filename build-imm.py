#! /usr/bin/env python

'''
This file identifies 
'''

import sys
import math
import scipy.stats
import numpy as np
import math
import compareResults
import time
import argparse

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

trueFile = ''

# If a kmer has at least this many occurences, don't use a smaller kmer
countThreshold = 400

def init():
    ''' Initialize counts and immScores tables to store training counts and calculated IMM probabilities.
    '''
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
    ''' Read a genome in FASTA format and return it as a string.
    '''
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def findLongORFs(genome):
    ''' Find long potential ORFs in the genome 
        (regions that start with a start codon and end with a stop codon in the same reading frame)
        Any of these potential ORFs that do not overlap are very likely to be true ORFs
    '''

    minLen = 500
    counts = [0,0]
    longORFs = []
    startPos = None

    # Find forward ORFs
    for i in xrange(3):
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue

            # Each ORF is defined by the first start codon and first stop codon that appear in the same reading frame
            if startPos == None and codon in start:
                startPos = nt
            elif not startPos == None and codon in stop:
                if nt+3-startPos >= minLen:
                    longORFs += [(startPos, nt+3, 0)]
                    #orfs += [(genome[startPos : nt+3], 0)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    # Find reverse complemented ORFs
    for i in xrange(3):
        startPos = None
        for pos in xrange(i, len(genome), 3):
            nt = len(genome)-pos
            codon = genome[nt-3:nt]
            if len(codon) < 3:
                continue

            # Each ORF is defined by the first start codon and first stop codon that appear in the same reading frame
            if startPos == None and codon in revCompStart:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                if startPos-nt+3 >= minLen:
                    longORFs += [(nt-3, startPos, 3)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                startPos = None

    # Remove any overlapping ORFs
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

        # Noncoding segments must start after a stop codon and end before a start codon
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
    ''' Update the lists of coding and noncoding regions based on the list of orfs passed in
    '''
    minLen = 500
    orfs = sorted(orfs)
    noncoding = []

    # noncoding regions must have some minimum length and lie between 2 orfs
    for i in xrange(len(orfs)-1):
        if orfs[i+1][0] - orfs[i][1] >= minLen:
            noncoding += [genome[orfs[i][1]:orfs[i+1][0]]]

    # all orfs passed in are coding sequences
    coding = []
    for c in orfs:
        coding += [(genome[c[0]:c[1]], c[2])]

    return coding, noncoding

def buildIMM(counts, codingSeqs, noncodingSeqs):
    ''' Trains a fixed-length Markov model

        counts: Empty list to be filled with kmer counts
        codingSeqs: ORFs for training
        noncodingSeqs: Regions that are definitely not ORFs, for training reading frame 6
    '''

    # Train reading frames 0-5 on  ORFs
    for (seq,frame) in codingSeqs:
        # frame = 0 (forward) or 3 (reverse complement)
        rf = frame
        for i in xrange(len(seq)):
            for l in xrange(min(maxLength+1, i+1)):
                kmer = seq[i-l:i+1]
                if kmer in counts[rf][l]:
                    counts[rf][l][kmer] += 1
                else:
                    counts[rf][l][kmer] = 1
            rf = (rf+1) % 3 + frame

    # Train reading frame 6 on non-ORFs
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
    ''' Trains a fixed-length Markov model

        counts: Empty list to be filled with kmer counts
        codingSeqs: ORFs for training
        noncodingSeqs: Regions that are definitely not ORFs, for training reading frame 6
        length: kmer context length
    '''

    # Train reading frames 0-5 on ORFs
    for (seq,frame) in codingSeqs:
        # frame = 0 (forward) or 3 (reverse complement)
        rf = (length % 3) + frame
        for i in xrange(length, len(seq)):
            kmer = seq[i-length:i+1]
            if kmer in counts[rf][length]:
                counts[rf][length][kmer] += 1
            else:
                counts[rf][length][kmer] = 1
            rf = (rf+1) % 3 + frame

    # Train reading frame 6 on non-ORFs
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
        Use the IMM or fixed-length Markov model to find all ORFs in the given genome

        genome: genome for which to find ORFs
        counts: list of kmer counts computed from training data
        immScore: list of immScores for kmers (filled in lazily to save time)
        mcLength: None for IMM, set to context length for fixed-length Markov Chain
    '''

    ''' Find potential ORFs in each '''
    orfs = []
    for i in xrange(6):
        orfs.append([])

    # Find forward ORFs
    for i in xrange(3):
        startPos = None
        for nt in xrange(i, len(genome), 3):
            codon = genome[nt:nt+3]
            if len(codon) < 3:
                continue
                
            # Each ORF is defined by the first start codon and first stop codon that appear in the same reading frame
            if codon in start and startPos == None:
                startPos = nt
            elif not startPos == None and codon in stop:
                if nt+3-startPos >= orfMinLength:
                    orfs[i] += [(startPos, nt+3)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                # Some nucleotides are unknown and appear as different letters; ignore them
                startPos = None

    # Find reverse complemented ORFs
    for i in xrange(3):
        startPos = None
        for pos in xrange(i, len(genome), 3):
            nt = len(genome)-pos
            codon = genome[nt-3:nt]
            if len(codon) < 3:
                continue

            # Each ORF is defined by the first start codon and first stop codon that appear in the same reading frame
            if codon in revCompStart and startPos == None:
                startPos = nt
            elif not startPos == None and codon in revCompStop:
                if startPos-nt+3 >= orfMinLength:
                    orfs[i+3] += [(nt-3, startPos)]
                startPos = None
            elif (codon[0] not in nts) or (codon[1] not in nts) or (codon[2] not in nts):
                # Some nucleotides are unknown and appear as different letters; ignore them
                startPos = None

    ''' filter out bad-scoring orfs '''
    goodORFs = []

    totals = [0]*6
    approved = [0]*6

    for i in xrange(len(orfs)):
        if debug:
            print 'Processing rf %d (%d)' % (i, len(orfs[i]))
        
        for orf in orfs[i]:
            # Score ORF in all reading frames
            if mcLength == None:
                rfScores = scoreIMM(orf, counts, immScores)
            else:
                rfScores = scoreMarkovChain(orf, mcLength, counts)

            # Only keep ORF if highest score is in the correct reading frame,
            #  and is greater than some threshold.
            maxScore = max(rfScores)
            length = orf[1]-orf[0]
            if i < 3 and maxScore == rfScores[0] and rfScores[0]/length > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[0]/length))
                approved[i] += 1
            elif i >= 3 and maxScore == rfScores[3] and rfScores[3]/length > scoreThreshold:
                goodORFs.append((orf[0], orf[1], i, rfScores[3]/length))
                approved[i] += 1
            
            totals[i] += 1

    ''' resolve overlapping orfs '''
    minOverlap = 20
    goodORFs = sorted(goodORFs)
    approvedORFs = []
    suspectORFs = []

    # If 2 ORFs overlap, score the overlapping region in each reading frame.
    # If one ORF is longer than the other and the highest scoring reading frame of the overlap
    #  is the same as the reading frame of the longer ORF, discard the shorter ORF.
    i = 0
    while i < len(goodORFs)-1:
        j = i+1
        while j < len(goodORFs) and (goodORFs[i][1] - goodORFs[j][0] >= minOverlap):
            if mcLength == None:
                overlapScores = scoreIMM((goodORFs[j][0], goodORFs[i][1]), counts, immScores)
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

    # Any remaining ORFs that overlap are suspect. All others are approved.
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
    ''' Write the set of ORFs to a file 
    '''

    print '  %d approved, %d suspected ORFs found' % (len(approvedORFs), len(suspectORFs))

    with open('predicted.txt', 'w') as f:
        for orf in approvedORFs:
            f.write(str(orf[0]) + '\t' + str(orf[1]) + '\n')
        for orf in suspectORFs:
            f.write('*' + str(orf[0]) + '\t' + str(orf[1]) + '\n')

def scoreIMM(orf, counts, immScores):
    ''' Return a list of the IMM scores for the given orf in each of the 7 reading frames.
        counts is the set of kmer counts computer from the training data.
    '''

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
    ''' Return a list of the fixed-length Markov model scores for the given orf in each of the 7 reading frames.
        counts is the set of kmer counts computer from the training data.
    '''
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
    ''' Return the interpolated Markov model score for the last character of seq appearing in the given
        reading frame following the first k-1 characters of seq.

        counts contains the counts of each kmer in each reading frame, calculated during training
        immScores contains the probabilities of each kmer in each reading frame (computed lazily)
    '''

    # Check if we have calculated this IMM score yet
    if seq in immScores[rf][len(seq)-1]:
        return immScores[rf][len(seq)-1][seq]

    prefix = seq[:-1]

    # If the sequence is of length 1 (i.e. no previous bases), we cannot interpolate any further, so IMM = Pr
    if len(seq) == 1:
        probs = [0]*4
        for i in xrange(4):
            immScores[rf][0][prefix+nts[i]] = prob(rf, prefix+nts[i], counts)
        return immScores[rf][len(seq)-1][seq]


    length = len(seq)-1
    if prefix in counts[rf][length-1] and counts[rf][length-1][prefix] > countThreshold:
        # If we see the context at least countThreshold times, we can be very confident in the probability, 
        #  so don't use any shorter contexts
        wgt = 1
    else:
        # Compute the weight to assign to this context by calculating how similar the set of counts for
        #  this context is to the IMM scores for the next smallest context
        currCounts = [1]*4
        currTotal = 4
        for i in xrange(4):
            tempSeq = prefix + nts[i]
            if tempSeq in counts[rf][length]:
                currCounts[i] += counts[rf][length][tempSeq]
                currTotal += counts[rf][length][tempSeq]

        nextProbs = [0]*4
        for i in xrange(4):
            nextProbs[i] = imm(rf, prefix[1:]+nts[i], counts, immScores)
        
        d = scipy.stats.chi2_contingency(np.array([currCounts, nextProbs]))[1]

        if d < 0.5:
            wgt = 0
        else:
            if prefix in counts[rf][length-1]:
                wgt = d * counts[rf][length-1][prefix] / 400.0
            else:
                wgt = 0

    # Compute the IMM score as a weight sum of the probability for this context + the IMM score for the next smallest context
    probs = [0]*4
    for i in xrange(4):
        probs[i] = wgt * prob(rf, prefix+nts[i], counts) + (1-wgt) * imm(rf, prefix[1:]+nts[i], counts, immScores)
    totalProb = sum(probs)

    # Normalize IMM scores to 1
    for i in xrange(4):
        immScores[rf][len(seq)-1][prefix+nts[i]] = float(probs[i]) / totalProb

    return immScores[rf][len(seq)-1][seq]

def prob(rf, seq, counts):
    ''' Return the probability of the last character of seq appearing in the given reading frame following 
        the first k-1 characters of seq. Applies +1 smoothing to the counts.
    '''
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
    ''' Run the IMM classifier and print the accuracy of the results.
    '''
    print 'IMM max length %d:' % maxLength
    counts, immScores = init()

    startTime = time.time()
    coding, noncoding = findLongORFs(genome)
    counts = buildIMM(counts, coding, noncoding)
    buildTime = time.time()
    approved, suspect = glimmer(genome, counts, immScores)
    endTime = time.time()
    print '  %0.2fs total (%0.2fs to build, %0.2fs to run)' % (endTime-startTime, buildTime-startTime, endTime-buildTime)
    writeORFs(approved, suspect)

    compareResults.compare(trueFile, 'predicted.txt')
    print ''

def runMC(length):
    ''' Run the fixed-length Markov model classifier and print the accuracy of the results.
    '''
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

    compareResults.compare(trueFile, 'predicted.txt')
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
        compareResults.compare(trueFile, 'predicted.txt')
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
    ''' Repeatedly train on the set of predicted ORFs, and use a fixed-length MM to predict a new set of ORFs.
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
        compareResults.compare(trueFile, 'predicted.txt')
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
    ''' Returns true if the sets of ORFs are exactly the same, false otherwise.
    '''
    if not len(oldORFs) == len(newORFs):
        return False

    oldORFs.sort()
    newORFs.sort()
    return (oldORFs == newORFs)

if __name__ == '__main__':
    ''' Main function. Runs either IMM or fixed-length MC classifier, iterative or non-iterative.
    '''


    # Print file's docstring if -h is invokedc
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genome', type=str, required=True, 
        help='Full path to FASTA file containing genome to process')
    parser.add_argument('--max-length', type=int, required=True, 
        help='Maximum kmer length')
    parser.add_argument('--trueORFs', type=str, required=True, 
        help='File containing true gene boundaries')
    parser.add_argument("--fixed", help="Run fixed-length Markov Chain. If not present, run IMM",
        action="store_true")
    parser.add_argument("--iterative", help="Run iterative. If not present, run once.",
        action="store_true")
    
    args = parser.parse_args(sys.argv[1:])
    
    genome = readFASTA(args.genome)
    maxLength = args.max_length
    trueFile = args.trueORFs
    
    if args.iterative:
        if args.fixed:
            runIterativeMC(maxLength)
        else:
            runIterativeIMM()

    else:
        if args.fixed:
            runMC(maxLength)
        else:
            runIMM()





    