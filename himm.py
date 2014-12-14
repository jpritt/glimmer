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
    
    countThreshold = 400.0

    def __init__(self, counts, singleProbs, probSwitch, maxLength):
        ''' Initialize the HIMM given 
            counts - counts of each kmer
            immProbs - empty, to be filled out with calculated IMM probs
            singleProbs - For each state, probability of each single nucleotide occurring
            maxLength - maximum kmer length used
        '''
        
        self.counts = counts
        self.singleProbs = singleProbs
        self.probSwitch = probSwitch
        self.maxLength = maxLength

        self.nts = ['A', 'C', 'G', 'T']
        self.stateNTs = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']

        # For each state, possible next states with non-zero probability
        self.nextStates = [[1,3], [2,3], [0,3], [1,3], [2,3], [0,3]]
        #self.nextStates = [[1,3], [2,4], [0,5], [1,3], [2,4], [0,5]]

        '''
        self.immProbs = []
        for i in xrange(len(counts)):
            self.immProbs.append([])
            for j in xrange(len(counts[0])):
                self.immProbs[i].append([])
                for k in xrange(self.maxLength+1):
                    self.immProbs[i][j].append(dict())
        '''

        # 6 tables: 0 -> (1,3), 1 -> (2,4), 2 -> (0,5), 3 -> (1,3), 4 -> (2,4), 5 -> (0,5)
        self.immProbs = []
        for i in xrange(6):
            self.immProbs.append([])
            for j in xrange(self.maxLength+1):
                self.immProbs[i].append(dict())
    
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
                    probA = mat[i-1, (newState-1)%3] + math.log(self.imm((newState-1)%3, newState, x[j-length:j+1]), 2) + self.probSwitch[0][0]
                    probB = mat[i-1, (newState-1)%3+3] + math.log(self.imm(3, newState, x[j-length:j+1]), 2) + self.probSwitch[1][0]
                    if probA > probB:
                        mat[i, j] = probA
                        matTb[i, j] = (newState-1)%3
                    else:
                        mat[i, j] = probB
                        matTb[i, j] = (newState-1)%3 + 3

                # Intronic states
                else:
                    # possible previous states: newState-3, newState
                    probA = mat[i-1, newState-3] + math.log(self.imm(newState-3, 3, x[j-length:j+1]), 2) + self.probSwitch[0][1]
                    probB = mat[i-1, newState] + math.log(self.imm(3, 3, x[j-length:j+1]), 2) + self.probSwitch[0][0]
                    if probA > probB:
                        mat[i, j] = probA
                        matTb[i, j] = newState-3
                    else:
                        mat[i, j] = probB
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
        if newState > 2:
            seq = seq[:-1] + seq[-1].lower()
        if seq in self.immProbs[oldState][len(seq)-1]:
            return self.immProbs[oldState][len(seq)-1][seq]

        prefix = seq[:-1]

        if len(seq) == 1:
            probs = [0]*4
            for i in xrange(4):
                self.immProbs[oldState][0][self.nts[i]] = self.singleProb(newState, self.nts[i])
            return self.immProbs[oldState][0][seq]

        length = len(seq)-1
        if prefix in self.counts[oldState][length-1] and self.counts[oldState][length-1][prefix] > self.countThreshold:
            wgt = 1
        else:
            probs = self.getProbs(oldState, prefix)

            currProbs = []
            nextProbs = []
            for n in self.stateNTs:
                currProbs.append(probs[n])
                nextProbs.append(self.imm(oldState, newState, prefix[1:]+n))
            
            d = scipy.stats.chi2_contingency(np.array([currCounts, nextProbs]))[1]

            if d < 0.5:
                wgt = 0
            else:
                if prefix in self.counts[oldState][newState][length-1]:
                    wgt = d * self.counts[oldState][newState][length-1][prefix] / self.countThreshold
                else:
                    wgt = 0

        probs = self.getProbs(oldState, prefix)
        for k in probs.keys():
            probs[i] = wgt * probs[k] + (1-wgt) * self.imm(oldState, newState, prefix[1:]+self.nts[i])
        totalProb = sum(probs)

        for i in xrange(4):
            self.immProbs[oldState][newState][len(seq)-1][prefix+self.nts[i]] = float(probs[i]) / totalProb

        return self.immProbs[oldState][newState][len(seq)-1][seq]

    def getProbs(self, oldState, prefix):
        #if newState > 3:
        #    newState = 3

        #if not newState in self.nextStates[oldState]:
        #    return None

        length = len(seq)-1

        # Count all possible next nucleotides, with +1 smoothing
        counts = dict()
        totalCount = 0
        for i in self.nextStates[oldState]:
            for n in self.nts:
                if i < 3:
                    stateN = n
                else:
                    stateN = n.lower()
                
                if seq[:-1]+n in self.counts[oldState][i][length]:
                    counts[stateN] =  1 + self.counts[oldState][i][length][seq[:-1]+n]
                else:
                    counts[stateN] =  1
                totalCount += counts[stateN]

        for k in counts.keys():
            counts[k] = float(counts[k]) / totalCount
        return probs

    '''
    def prob(self, oldState, newState, seq):
        if newState > 3:
            newState = 3

        if not newState in self.nextStates[oldState]:
            return 0

        length = len(seq)-1

        # Count all possible next nucleotides, with +1 smoothing
        totalNum = len(self.stateNTs)
        for i in self.nextStates[oldState]:
            for n in self.nts:
                if seq[:-1]+n in self.counts[oldState][newState][length]:
                    totalNum += self.counts[oldState][newState][length][seq[:-1]+n]

        if not seq in self.counts[oldState][newState][length]:
            return 1.0 / float(totalNum)
        else:
            return float(1+self.counts[oldState][newState][length][seq]) / float(totalNum)
    '''

    def singleProb(self, state, nt):
        totalNum = 4
        for n in self.nts:
            totalNum += self.singleProbs[state][n]

        return float(1+self.singleProbs[state][nt]) / float(totalNum)