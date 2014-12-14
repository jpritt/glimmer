#! /usr/bin/env python

class HIMMTrain:

    def __init__(maxLen):
        self.maxLength = maxLen

        # 4 states: intron, exon_rf1, exon_rf2, exon_rf3
        # counts[i][j] contains a dictionary with kmer counts from state i to state j
        self.counts = []
        for i in xrange(4):
            self.counts.append([])
            for j in xrange(4):
                self.counts[i].append([])
                for k in xrange(self.maxLength+1):
                    self.counts[i][j].append(dict())

    def train(gene, exonMask):
        '''
            gene - string of nucleotides (A,C,G,T only)
            exonMask - string of 1s (exon) and 0s (intron) of same length as genome
        '''

