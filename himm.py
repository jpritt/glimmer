#! /usr/bin/env python
import math
import numpy as np
import scipy.stats
import sys

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

        self.interpolated = True
        self.fixedLength = 0

    def setFixedLength(self, length):
        self.interpolated = False
        self.fixedLength = length

    def setInterpolated(self):
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

            #if j > 885 and j < 890:
            #    print 'j = %d' % j

            for newState in xrange(0, nrow):
                # Exonic states
                if newState < 3:
                    # possible previous states: (newState-1)%3, (newState-1)%3 + 3
                    #if self.interpolated:
                    #print 'Calculating imm(%d)' % (x[j-length:j+1])
                    probA = mat[(newState-1)%3, j-1] + math.log(self.imm((newState-1)%3, newState, x[j-length:j+1]), 2)
                    probB = mat[(newState-1)%3+3, j-1] + math.log(self.imm(3, newState, x[j-length:j+1]), 2)
                    #else:
                    #    probsA = self.getProbs((newState-1)%3, x[j-length:j])
                    #    probA = mat[(newState-1)%3, j-1] + math.log(probsA[x[j]], 2)

                    #    probsB = self.getProbs(3, x[j-length:j])
                    #    probB = mat[(newState-1)%3+3, j-1] + math.log(probsB[x[j]], 2)

                    #print '  State %d:' % newState
                    #print '    %d --> %f (%f + %f)' % ((newState-1)%3, probA, mat[(newState-1)%3, j-1], math.log(self.imm((newState-1)%3, newState, x[j-length:j+1]), 2))
                    #print '    %d --> %f (%f + %f)' % ((newState-1)%3+3, probB, mat[(newState-1)%3+3, j-1], math.log(self.imm(3, newState, x[j-length:j+1]), 2))
                    if probA > probB:
                        mat[newState, j] = probA
                        matTb[newState, j] = (newState-1)%3
                    else:
                        mat[newState, j] = probB
                        matTb[newState, j] = (newState-1)%3 + 3

                # Intronic states
                else:
                    # possible previous states: newState-3, newState
                    #if self.interpolated:
                    probA = mat[newState-3, j-1] + math.log(self.imm(newState-3, 3, x[j-length:j+1]), 2)
                    probB = mat[newState, j-1] + math.log(self.imm(3, 3, x[j-length:j+1]), 2)
                    #else:
                    #    probsA = self.getProbs(newState-3, x[j-length:j])
                    #    probA = mat[newState-3, j-1] + math.log(probsA[x[j].lower()], 2)

                    #    probsB = self.getProbs(3, x[j-length:j])
                    #    probB = mat[newState, j-1] + math.log(probsB[x[j].lower()], 2)


                    #print '  State %d:' % newState
                    #print '    %d --> %f' % (newState-3, probA)
                    #print '    %d --> %f' % (newState, probB)

                    if probA > probB:
                        mat[newState, j] = probA
                        matTb[newState, j] = newState-3
                    else:
                        mat[newState, j] = probB
                        matTb[newState, j] = newState

        '''
        for i in xrange(6):
            for j in xrange(16,26):
                sys.stdout.write('%0.2f\t' % mat[i,j])
                #sys.stdout.write(str(mat[i,j]) + '\t')
            print ''
        print ''
        for i in xrange(6):
            for j in xrange(16,26):
                sys.stdout.write('%d\t' % matTb[i,j])
                #sys.stdout.write(str(matTb[i,j]) + '\t')
            print ''
        print ''
        exit()
        '''
        

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
            probs = self.getProbs(oldState, '')
            totalProb = 0
            for k in probs.keys():
                totalProb += probs[k]

            for n in self.stateNTs:
                self.immProbs[oldState][0][n] = float(probs[n]) / totalProb

            return self.immProbs[oldState][0][seq]

        length = len(seq)-1

        if self.interpolated:
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
            wgt = 1

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

        for n in self.stateNTs:
            self.immProbs[oldState][len(seq)-1][prefix+n] = float(probs[n]) / totalProb

        return self.immProbs[oldState][len(seq)-1][seq]

    def getProbs(self, oldState, prefix):
        #if newState > 3:
        #    newState = 3

        #if not newState in self.nextStates[oldState]:
        #    return None

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