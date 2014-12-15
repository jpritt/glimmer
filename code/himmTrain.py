#! /usr/bin/env python
import math

class HIMMTrain:
    ''' Trains a HIMM model.
    '''

    def __init__(self, maxLen):
        self.maxLength = maxLen

        # Initialize counts tables
        # 4 states: intron, exon_rf1, exon_rf2, exon_rf3
        # counts[i][j] contains a dictionary with kmer counts from state i to state j
        self.counts = []
        for i in xrange(4):
            self.counts.append([])
            for j in xrange(4):
                self.counts[i].append([])
                for k in xrange(self.maxLength+1):
                    self.counts[i][j].append(dict())

    def train(self, gene, exonMask):
        '''
            Train the lists of counts on a single gene.

            gene: string of nucleotides (A,C,G,T only)
            exonMask: string of 1s (exon) and 0s (intron) of same length as genome
        '''

        rf = 0
        for i in xrange(1, len(gene)):
            # Find previous state gene was in
            if exonMask[i-1] == '1':
                prevState = rf
            else:
                prevState = 3

            # Find current state gene is in
            if exonMask[i] == '1':
                rf = (rf+1) % 3
                newState = rf
            else:
                newState = 3

            # Update appropriate counts table for each kmer ending at the current base
            for l in xrange(min(i, self.maxLength+1)):
                seq = gene[i-l:i+1]
                if seq in self.counts[prevState][newState][l]:
                    self.counts[prevState][newState][l][seq] += 1
                else:
                    self.counts[prevState][newState][l][seq] = 1

