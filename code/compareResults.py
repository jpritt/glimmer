#! /usr/bin/env python
''' Contains methods for comparing 2 sets of ORFs to calculate prediciton accuracy.
'''

import sys

def matches(orfA, orfB):
    ''' Since starting sites are hard to localize, 2 ORFs are considered the same if they begin or end on the same base.
        Salzberg et al. in the GLIMMER paper used the same method for measuring accuracy.
    '''
    if orfA[0] == orfB[0] or orfA[1] == orfB[1]:
        return True
    else:
        return False

def compare(trueORFs, predORFs):
    ''' Compare the 2 lists of ORFs and print the number of true positives, false positives, and false negatives.
    '''
    true = []
    with open(trueORFs, 'r') as f:
        for line in f:
            orf = line.rstrip().split('\t')
            true.append((int(orf[0]), int(orf[1])))

    pred = []
    susp = []
    with open(predORFs, 'r') as f:
        for line in f:
            if line[0] == '*':
                orf = line[1:].rstrip().split('\t')
                susp.append((int(orf[0]), int(orf[1])))
            else:
                orf = line.rstrip().split('\t')
                pred.append((int(orf[0]), int(orf[1])))


    # Count true positives, false positives, and false negatives    
    pred = pred + susp
    tp = 0
    fp = 0
    fn = 0
    for orf in pred:
        found = False
        for o in true:
            if matches(orf, o):
                tp += 1
                found = True
                break
        if not found:
            fp += 1
    for orf in true:
        found = False
        for o in pred:
            if matches(orf, o):
                found = True
                break
        if not found:
            fn += 1

    print '    %% found: %0.2f (%d)' % (float(100*tp)/len(true), tp)
    print '    Genes missed: %d' % fn
    print '    Additional genes: %d' % fp
    

