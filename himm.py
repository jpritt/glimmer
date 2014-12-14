#! /usr/bin/env python
import math
import numpy

class HIMM(object):
    ''' Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities in dictionaries
        mapping 2-character codes onto floating-point probabilities
        for those table entries.  States and emissions are represented
        with single characters.  Emission symbols comes from a finite.  '''
    
    def __init__(self, counts, immProbs, maxLength):
        ''' Initialize the HIMM given 
            E - Emission probabilities
            I - Initial probabilities
            counts - counts of each kmer
        '''
        
        self.counts = counts
        self.immProbs = immProbs
        self.maxLength = maxLength
    
    def jointProbL(self, p, x):
        ''' Return log2 of joint probability of path p and emission
            string x.'''
        p = map(self.qmap.get, p) # turn state characters into ids
        x = map(self.smap.get, x) # turn emission characters into ids
        tot = self.Ilog[p[0]] # start with initial probability
        for i in xrange(1, len(p)):
            tot += self.Alog[p[i-1], p[i]] # transition probability
        for i in xrange(0, len(p)):
            tot += self.Elog[p[i], x[i]] # emission probability
        return tot
    
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