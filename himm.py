#! /usr/bin/env python
import math
import numpy as np
import scipy.stats
import sys

class HIMM(object):
    ''' Interpolated Hidden Markov Model implementation.
    '''
    
    countThreshold = 400.0

    def __init__(self, counts, maxLength):
        ''' Initialize the HIMM. 

            counts - counts of kmers for each state transition
            maxLength - maximum kmer length used
        '''
        
        self.counts = counts
        self.maxLength = maxLength

        self.nts = ['A', 'C', 'G', 'T']

        # Each state can be follow by 1 of the exon states or 1 of the intron states.
        # Uppercase nts indicate a base in the exon state, lowercase nts indicate a base in the intron state.
        self.stateNTs = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']

        # For each state, possible next states with non-zero probability
        self.nextStates = [[1,3], [2,3], [0,3], [1,3], [2,3], [0,3]]

        # Tables containing IMM probabilities, one for each possible state transiton.
        # Tables are initially empty, filled lazily to save time.
        # 6 tables: 0 -> (1,3), 1 -> (2,4), 2 -> (0,5), 3 -> (1,3), 4 -> (2,4), 5 -> (0,5)
        self.immProbs = []
        for i in xrange(6):
            self.immProbs.append([])
            for j in xrange(self.maxLength+1):
                self.immProbs[i].append(dict())

        # If not interpolated, uses a standard fixed-length Hidden Markov Model
        self.interpolated = True
        self.fixedLength = 0

    def setFixedLength(self, length):
        ''' Use a fixed-length HMM. '''
        self.interpolated = False
        self.fixedLength = length

    def setInterpolated(self):
        ''' Use an Interpolated HMM. '''
        self.interpolated = True
    
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
        if self.interpolated:
            mat[0, 0] = 0
            for i in xrange(1, nrow):
                mat[i, 0] = -100000
            startBase = 1
        else:
            # If fixed-length, set first fixedLength states to be exons.
            for j in xrange(self.fixedLength):
                mat[0, j] = 0
                for i in xrange(1, nrow):
                    mat[i, j] = -100000
            startBase = self.fixedLength

        # Fill in rest of log prob and Tb tables
        for j in xrange(startBase, ncol):
            if self.interpolated:
                length = min(j, self.maxLength)
            else:
                length = self.fixedLength

            for newState in xrange(0, nrow):
                # Exonic states
                if newState < 3:
                    probA = mat[(newState-1)%3, j-1] + math.log(self.imm((newState-1)%3, newState, x[j-length:j+1]), 2)
                    probB = mat[(newState-1)%3+3, j-1] + math.log(self.imm(3, newState, x[j-length:j+1]), 2)
                    
                    if probA > probB:
                        mat[newState, j] = probA
                        matTb[newState, j] = (newState-1)%3
                    else:
                        mat[newState, j] = probB
                        matTb[newState, j] = (newState-1)%3 + 3

                # Intronic states
                else:
                    probA = mat[newState-3, j-1] + math.log(self.imm(newState-3, 3, x[j-length:j+1]), 2)
                    probB = mat[newState, j-1] + math.log(self.imm(3, 3, x[j-length:j+1]), 2)

                    if probA > probB:
                        mat[newState, j] = probA
                        matTb[newState, j] = newState-3
                    else:
                        mat[newState, j] = probB
                        matTb[newState, j] = newState
        

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
        ''' Return the interpolated Markov model score for the last character of seq appearing in newState
            following the first k-1 characters of seq in oldState.
        '''

        # Use lowercase characters for intron regions
        if newState > 2:
            seq = seq[:-1] + seq[-1].lower()

        # If it has already been calculated, no need to recomputer
        if seq in self.immProbs[oldState][len(seq)-1]:
            return self.immProbs[oldState][len(seq)-1][seq]


        prefix = seq[:-1]

        # If length is 1, don't interpolate. Just calculate P for each character that could come after oldState.
        if len(seq) == 1:
            probs = self.getProbs(oldState, '')
            totalProb = 0
            for k in probs.keys():
                totalProb += probs[k]

            for n in self.stateNTs:
                self.immProbs[oldState][0][n] = float(probs[n]) / totalProb

            return self.immProbs[oldState][0][seq]

        length = len(seq)-1
        if self.interpolated:
            # Assign a weight to the probability for the current context length based on the Chi-square
            #  value of the probabilities for the current context length compared to the interpolated
            #   probabilities for the next smallest context length.
            relCount = 0
            for s in self.nextStates[oldState]:
                if prefix in self.counts[oldState][s][len(prefix)-1]:
                    relCount += self.counts[oldState][s][len(prefix)-1][prefix]

            if relCount > self.countThreshold:
                wgt = 1
            else:
                probs = self.getProbs(oldState, prefix)

                currProbs = []
                nextProbs = []
                for n in self.stateNTs:
                    currProbs.append(probs[n])

                    if n.islower():
                        newS = self.nextStates[oldState][1]
                    else:
                        newS = self.nextStates[oldState][0]
                    nextProbs.append(self.imm(oldState, newS, prefix[1:]+n))
                
                d = scipy.stats.chi2_contingency(np.array([currProbs, nextProbs]))[1]

                if d < 0.5:
                    wgt = 0
                else:
                    wgt = d * relCount / self.countThreshold
        else:
            # If we're using the fixed-length model, don't interpolate at all, just return P
            wgt = 1

        # A character in oldState can be followed by 'A', 'C', 'G', or 'T' in a certain exonic newState, or
        #   'a', 'c', 'g', or 't', in a certain intronic newState.
        # Compute the interpolated probabilities of each of these and store in the immProbs table for oldState.
        probs = self.getProbs(oldState, prefix)
        totalProb = 0
        for k in probs.keys():
            if k.islower():
                newS = self.nextStates[oldState][1]
            else:
                newS = self.nextStates[oldState][0]

            if wgt == 0:
                probs[k] = self.imm(oldState, newS, prefix[1:]+k)
            elif wgt < 1:
                probs[k] = wgt * probs[k] + (1-wgt) * self.imm(oldState, newS, prefix[1:]+k)
            totalProb += probs[k]

        # Normalize probabilities to 1
        for n in self.stateNTs:
            self.immProbs[oldState][len(seq)-1][prefix+n] = float(probs[n]) / totalProb

        return self.immProbs[oldState][len(seq)-1][seq]

    def getProbs(self, oldState, prefix):
        ''' Return the list of probabilities of each nucleotide appearing after the given prefix in the given oldState.
            Uppercase bases indicate that the new state is the exonic reading frame state which can appear after oldState.
            Lowercase bases indicate that the new state is the intron reading frame state which can appear after oldState.
        '''

        length = len(prefix)

        # Count all possible next nucleotides, with +1 smoothing
        counts = dict()
        totalCount = 0
        for i in self.nextStates[oldState]:
            for n in self.nts:
                if i < 3:
                    stateN = n
                else:
                    stateN = n.lower()

                if (prefix+n) in self.counts[oldState][i][length]:
                    counts[stateN] =  1 + self.counts[oldState][i][length][prefix+n]
                else:
                    counts[stateN] =  1
                totalCount += counts[stateN]

        for k in counts.keys():
            counts[k] = float(counts[k]) / totalCount
        return counts