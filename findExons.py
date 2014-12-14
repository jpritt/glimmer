#! /usr/bin/env python
import sys
import himm
import himmTrain

def runHIMM(genome, trainGenes, testGenes, maxLength):
    trainer = himmTrain.HIMMTrain(maxLength)

    for gene in trainGenes:
        trainer.train(genome[gene[0]:gene[1]], gene[2])
    print trainer.singleProbs

    model = himm.HIMM(trainer.counts, trainer.singleProbs, maxLength)

    for gene in testGenes:
        predictedMask = model.viterbi(genome[gene[0]:gene[1]])

        print '%d, %d' % (len(predictedMask), len(gene[2]))

        correct = 0
        for i in xrange(len(predictedMask)):
            if predictedMask[i] == gene[2][i]:
                correct += 1

        print '%0.2f correct (%d / %d)' % (float(correct)/float(len(predictedMask)), correct, len(predictedMask))
        exit()

def readGenes(filename):
    genes = []
    with open(filename, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            bounds = []
            for r in row:
                b = r.split(',')
                bounds.append((int(b[0]), int(b[1])))

            start = bounds[0][0]
            end = bounds[-1][-1]
            mask = ''
            for i in xrange(len(bounds)):
                if i > 0:
                    mask += '0' * (bounds[i][0] - bounds[i-1][1])
                mask += '1' * (bounds[i][1] - bounds[i][0])

            genes.append((start, end, mask))
    return genes

def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

genome = readFASTA(sys.argv[1])
genes = readGenes(sys.argv[2])

trainSize = len(genes) / 10
runHIMM(genome, genes[:trainSize], genes[trainSize:], 6)