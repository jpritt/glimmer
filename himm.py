#! /usr/bin/env python
import math
import numpy

class HIMM(object):
    ''' Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities in dictionaries
        mapping 2-character codes onto floating-point probabilities
        for those table entries.  States and emissions are represented
        with single characters.  Emission symbols comes from a finite.  '''
    
    def __init__(self, counts, immProbs, singleProbs, maxLength):
        ''' Initialize the HIMM given 
            counts - counts of each kmer
            immProbs - empty, to be filled out with calculated IMM probs
            singleProbs - For each state, probability of each single nucleotide occurring
            maxLength - maximum kmer length used
        '''
        
        self.counts = counts
        self.immProbs = immProbs
        self.singleProbs = singleProbs
        self.maxLength = maxLength

        self.nts = ['A', 'C', 'G', 'T']
    
    def viterbiL(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with log2 of its probability.'''

        # 6 states: exon_rf1, exon_rf2, exon_rf3, intron_rf1, intron_rf2, intron_rf3
        nrow = 6
        ncol = len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # backtrace

        # Fill in first column
        # Every gene starts in exon_rf1
        mat[0, 0] = 1
        for i in xrange(1, nrow):
            mat[i, 0] = self.Elog[i, x[0]] + self.Ilog[i]

        # Fill in rest of log prob and Tb tables
        for j in xrange(1, ncol):
            for i in xrange(0, nrow):
                ep = self.Elog[i, x[j]]
                mx, mxi = mat[0, j-1] + self.Alog[0, i] + ep, 0
                for i2 in xrange(1, nrow):
                    pr = mat[i2, j-1] + self.Alog[i2, i] + ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi

        # Find final state with maximal log probability
        omx, omxi = mat[0, ncol-1], 0
        for i in xrange(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Backtrace
        i, p = omxi, [omxi]
        for j in xrange(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = ''.join(map(lambda x: self.Q[x], p[::-1]))
        return omx, p # Return log probability and path

    def imm(oldState, newState, seq):
        if seq in self.immProbs[oldState][newState][len(seq)-1]:
            return self.immProbs[oldState][newState][len(seq)-1]

        prefix = seq[:-1]

        if len(seq) == 1:
            probs = [0]*4
            for i in xrange(4):
                self.immProbs[oldState][newState][0][nts[i]] = self.singleProbs[newState][seq]
            return self.immProbs[oldState][newState][0][seq]

        length = len(seq)-1
        if prefix in self.counts[oldState][newState][length-1] and self.counts[oldState][newState][length-1][prefix] > countThreshold:
            wgt = 1
        else:
            currCounts = [1]*4
            currTotal = 4
            for i in xrange(4):
                tempSeq = prefix + nts[i]
                if tempSeq in self.counts[oldState][newState][length]:
                    currCounts[i] += self.counts[oldState][newState][length][tempSeq]
                    currTotal += self.counts[oldState][newState][length][tempSeq]

            nextProbs = [0]*4
            for i in xrange(4):
                nextProbs[i] = imm(oldState, newState, prefix[1:]+nts[i])
            
            d = scipy.stats.chi2_contingency(np.array([currCounts, nextProbs]))[1]

            if d < 0.5:
                wgt = 0
            else:
                if prefix in self.counts[oldState][newState][length-1]:
                    wgt = d * self.counts[oldState][newState][length-1][prefix] / 400.0
                else:
                    wgt = 0

        probs = [0]*4
        for i in xrange(4):
            probs[i] = wgt * prob(oldState, newState, prefix+nts[i]) + (1-wgt) * imm(oldState, newState, prefix[1:]+nts[i])
        totalProb = sum(probs)

        for i in xrange(4):
            self.immProbs[oldState][newState][len(seq)-1][prefix+nts[i]] = float(probs[i]) / totalProb

        return self.immProbs[oldState][newState][len(seq)-1][seq]

    def prob(oldState, newState, seq):
        length = len(seq)-1

        # Count all possible next nucleotides, with +1 smoothing
        totalNum = 0
        for n in nts:
            if seq[:-1]+n in self.counts[oldState][newState][length]:
                totalNum += self.counts[oldState][newState][length][seq[:-1]+n]
            else:
                totalNum += 1
        if not seq in self.counts[oldState][newState][length]:
            return 1.0 / float(totalNum)
        else:
            return float(self.counts[oldState][newState][length][seq]) / float(totalNum)