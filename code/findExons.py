#! /usr/bin/env python
import sys
import himm
import himmTrain
import time
import argparse

def runHIMM(genome, trainGenes, testGenes, maxLength, fixed):
    ''' Train a Hidden Interpolated Markov Model and use it to predict introns/exons for a set of genes.
    '''

    # Create trainer
    trainer = himmTrain.HIMMTrain(maxLength)

    # Train on training set
    for gene in trainGenes:
        trainer.train(genome[gene[0]:gene[1]], gene[2])

    # Create HIMM with trained data
    model = himm.HIMM(trainer.counts, maxLength)
    if fixed:
        model.setFixedLength(maxLength)

    # Test on testing set
    totalCorrect = 0
    totalBases = 0
    totalTime = 0
    totalLen = 0

    n = 0
    for gene in testGenes:
        # Predict a string of states
        startTime = time.time()
        predictedMask = model.viterbi(genome[gene[0]:gene[1]])
        endTime = time.time()

        totalTime += endTime - startTime
        totalLen += gene[1] - gene[0]

        # Compute accuracy of prediction as simple proportion of bases correct
        totalBases += len(gene[2])
        correct = 0
        for i in xrange(len(gene[2])):
            if predictedMask[i] == gene[2][i]:
                totalCorrect += 1
                correct += 1
        n += 1
        
        # Currently set to end after testing 1000 genes, for the same of time
        if n % 100 == 0:
            break

    print '  Final Accuracy: %0.2f' % (float(totalCorrect)/float(totalBases))
    print '  Average time per gene: %0.2f s' % (totalTime / n)
    print '  Average gene length: %0.1f' % (totalLen / n)

def readGenes(filename):
    ''' Reads gene locations from the given file.
    '''
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
    ''' Read a genome in FASTA format and return it as a string.
    '''
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

if __name__ == '__main__':
    ''' Main function. Runs either Hidden Interpolated Markov Model or fixed-length Hidden Markov Model.
        NOTE: Takes about a half hour to run.
    '''

    # Print file's docstring if -h is invokedc
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genome', type=str, required=True, 
        help='Path to FASTA file containing genome to process')
    parser.add_argument('--genes', type=str, required=True, 
        help='Path to file containing gene location information')
    parser.add_argument('--max-length', type=int, required=True, 
        help='Maximum kmer length')
    parser.add_argument("--fixed", help="Run fixed-length Markov Chain. If not present, run IMM",
        action="store_true")
    
    args = parser.parse_args(sys.argv[1:])
    
    genome = readFASTA(args.genome)
    genes = readGenes(args.genes)
    maxLength = args.max_length

    # Make training set first tenth of genes list
    trainSize = len(genes) / 10
    print 'Length = %d' % maxLength
    if args.fixed:
        print 'Fixed-length HMM:'
        runHIMM(genome, genes[:trainSize], genes[trainSize:], maxLength, True)
    else:
        print 'HIMM:'
        runHIMM(genome, genes[:trainSize], genes[trainSize:], maxLength, False)