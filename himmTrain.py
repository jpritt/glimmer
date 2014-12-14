#! /usr/bin/env python
import math

class HIMMTrain:

    def __init__(self, maxLen):
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

        self.singleProbs = []
        for i in xrange(4):
            self.singleProbs.append({'A':0, 'C':0, 'G':0, 'T':0})

        self.totalLen = [0,0]
        self.num = [0,0]

    def train(self, gene, exonMask):
        '''
            gene - string of nucleotides (A,C,G,T only)
            exonMask - string of 1s (exon) and 0s (intron) of same length as genome
        '''

        rf = 0
        start = 0
        for i in xrange(1, len(gene)):
            if not exonMask[i] == exonMask[i-1]:
                self.totalLen[int(exonMask[i-1])] += i - start
                self.num[int(exonMask[i-1])] += 1
                start = i

            if exonMask[i-1] == '1':
                prevState = rf
            else:
                prevState = 3

            if exonMask[i] == '1':
                rf = (rf+1) % 3
                newState = rf
            else:
                newState = 3

            for l in xrange(1, min(i, self.maxLength)):
                seq = gene[i-l:i+1]
                if seq in self.counts[prevState][newState][l]:
                    self.counts[prevState][newState][l][seq] += 1
                else:
                    self.counts[prevState][newState][l][seq] = 1

            self.singleProbs[newState][gene[i]] += 1

        self.totalLen[int(exonMask[-1])] += len(gene) - start
        self.num[int(exonMask[-1])] += 1

    def getProbSwitch(self):
        eToI = float(self.num[0]) / float(self.totalLen[0])
        iToE = float(self.num[1]) / float(self.totalLen[1])
        return [[math.log(1-eToI, 2), math.log(eToI, 2)], [math.log(iToE, 2), math.log(1-iToE, 2)]]
