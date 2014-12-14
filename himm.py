#! /usr/bin/env python
import math
import numpy as np
import scipy.stats

class HIMM(object):
    ''' Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities in dictionaries
        mapping 2-character codes onto floating-point probabilities
        for those table entries.  States and emissions are represented
        with single characters.  Emission symbols comes from a finite.  '''
    
    def __init__(self, counts, singleProbs, maxLength):
        ''' Initialize the HIMM given 
            counts - counts of each kmer
            immProbs - empty, to be filled out with calculated IMM probs
            singleProbs - For each state, probability of each single nucleotide occurring
            maxLength - maximum kmer length used
        '''
        
        self.counts = counts
        self.singleProbs = singleProbs
        self.maxLength = maxLength

        self.nts = ['A', 'C', 'G', 'T']

        self.immProbs = []
        for i in xrange(len(counts)):
            self.immProbs.append([])
            for j in xrange(len(counts[0])):
                self.immProbs[i].append([])
                for k in xrange(self.maxLength+1):
                    self.immProbs[i][j].append(dict())
    
    def viterbi(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with log2 of its probability.'''

        # 6 states: exon_rf1, exon_rf2, exon_rf3, intron_rf1, intron_rf2, intron_rf3
        nrow = 6
        ncol = len(x)
        mat   = np.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = np.zeros(shape=(nrow, ncol), dtype=int)   # backtrace

        # Fill in first column
        # Every gene starts in exon_rf1
        mat[0, 0] = 0
        for i in xrange(1, nrow):
            mat[i, 0] = -100000

        # Fill in rest of log prob and Tb tables
        for j in xrange(1, ncol):
            length = min(j, self.maxLength)
            for newState in xrange(0, nrow):
                # Exonic states
                if newState < 3:
                    # possible previous states: (newState-1)%3, (newState-1)%3 + 3
                    probA = mat[i-1, (newState-1)%3] + self.imm((newState-1)%3, newState, x[j-length:j+1])
                    probB = mat[i-1, (newState-1)%3+3] + self.imm(3, newState, x[j-length:j+1])
                    if probA > probB:
                        mat[i, j] = math.log2(probA)
                        matTb[i, j] = (newState-1)%3
                    else:
                        mat[i, j] = math.log2(probB)
                        matTb[i, j] = (newState-1)%3 + 3

                # Intronic states
                else:
                    # possible previous states: newState-3, newState
                    probA = mat[i-1, newState-3] + self.imm(newState-3, 3, x[j-length:j+1])
                    probB = mat[i-1, newState] + self.imm(3, 3, x[j-length:j+1])
                    if probA > probB:
                        mat[i, j] = math.log2(probA)
                        matTb[i, j] = newState-3
                    else:
                        mat[i, j] = math.log2(probB)
                        matTb[i, j] = newState

        # Find final state with maximal log probability
        maxProb = mat[0, ncol-1]
        maxState = 0
        for i in xrange(1, nrow):
            if mat[i, ncol-1] > maxProb:
                maxProb = mat[i, ncol-1]
                maxState = i

        # Backtrace
        if maxState < 3:
            path = '1'
        else:
            path = '0'
        for j in xrange(ncol-1, 0, -1):
            i = matTb[i, j]
            if i < 3:
                path = '1' + path
            else:
                path = '0' + path
        # Return path
        return path

    def imm(self, oldState, newState, seq):
        if seq in self.immProbs[oldState][newState][len(seq)-1]:
            return self.immProbs[oldState][newState][len(seq)-1][seq]

        prefix = seq[:-1]

        if len(seq) == 1:
            probs = [0]*4
            for i in xrange(4):
                self.immProbs[oldState][newState][0][self.nts[i]] = self.singleProb(newState, self.nts[i])
            return self.immProbs[oldState][newState][0][seq]

        length = len(seq)-1
        if prefix in self.counts[oldState][newState][length-1] and self.counts[oldState][newState][length-1][prefix] > countThreshold:
            wgt = 1
        else:
            currCounts = [1]*4
            currTotal = 4
            for i in xrange(4):
                tempSeq = prefix + self.nts[i]
                if tempSeq in self.counts[oldState][newState][length]:
                    currCounts[i] += self.counts[oldState][newState][length][tempSeq]
                    currTotal += self.counts[oldState][newState][length][tempSeq]

            nextProbs = [0]*4
            for i in xrange(4):
                nextProbs[i] = self.imm(oldState, newState, prefix[1:]+self.nts[i])
            
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
            probs[i] = wgt * self.prob(oldState, newState, prefix+self.nts[i]) + (1-wgt) * self.imm(oldState, newState, prefix[1:]+self.nts[i])
        totalProb = sum(probs)

        for i in xrange(4):
            self.immProbs[oldState][newState][len(seq)-1][prefix+self.nts[i]] = float(probs[i]) / totalProb

        return self.immProbs[oldState][newState][len(seq)-1][seq]

    def prob(self, oldState, newState, seq):
        length = len(seq)-1

        # Count all possible next nucleotides, with +1 smoothing
        totalNum = 4
        for n in self.nts:
            if seq[:-1]+n in self.counts[oldState][newState][length]:
                totalNum += self.counts[oldState][newState][length][seq[:-1]+n]

        if not seq in self.counts[oldState][newState][length]:
            return 1.0 / float(totalNum)
        else:
            return float(1+self.counts[oldState][newState][length][seq]) / float(totalNum)

    def singleProb(self, state, nt):
        totalNum = 4
        for n in self.nts:
            totalNum += self.singleProbs[newState][n]

        return float(1+self.singleProbs[newState][nt]) / float(totalNum)