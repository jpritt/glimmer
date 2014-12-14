#! /usr/bin/env python
import sys
import himm
import himmTrain
import time

def runHIMM(genome, trainGenes, testGenes, maxLength):
    trainer = himmTrain.HIMMTrain(maxLength)

    for gene in trainGenes:
        trainer.train(genome[gene[0]:gene[1]], gene[2])

    model = himm.HIMM(trainer.counts, trainer.singleProbs, trainer.getProbSwitch(), maxLength)

    totalCorrect = 0
    totalBases = 0
    totalTime = 0
    totalLen = 0

    n = 0
    for gene in testGenes:
        #print '%d - %d' % (n, gene[1]-gene[0])
        startTime = time.time()
        predictedMask = model.viterbi(genome[gene[0]:gene[1]])
        endTime = time.time()
        totalTime += endTime - startTime
        totalLen += gene[1] - gene[0]

        totalBases += len(gene[2])
        correct = 0
        for i in xrange(len(gene[2])):
            if predictedMask[i] == gene[2][i]:
                totalCorrect += 1
                correct += 1

        n += 1
        
        if True:#n % 10 == 0:
            #print 'Finished %d' % (n, len(testGenes))
            print '  Accuracy: %0.2f (%d, %0.2fs)' % (float(correct)/float(len(gene[2])), gene[1]-gene[0], endTime-startTime)
            #print '  Accuracy: %0.2f' % (float(totalCorrect)/float(totalBases))
            #print '  Time per gene: %0.2f s' % (totalTime / n)
            #print '  Average gene length: %0.1f' % (totalLen / n)
        if n >= 100:
            break

    print '\nFinal Accuracy: %0.2f' % (float(totalCorrect)/float(totalBases))
    print 'Average time per gene: %0.2f s' % (totalTime / len(testGenes))
    print 'Average gene length: %0.1f' % (totalLen / len(testGenes))

def readGenes(filename):
    genes = set()
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

            genes.add((start, end, mask))
    return list(genes)

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