#! /usr/bin/env python

import sys
def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def matches(orfA, orfB):
    if orfA[0] == orfB[0] or orfA[1] == orfB[1]:
        return True
    else:
        return False

    if orfA[0] < orfB[1] and orfA[1] > orfB[0]:
        return True
        olap = min(orfA[1],orfB[1]) - max(orfA[0],orfB[0])
        if olap > 0.75*(orfA[1]-orfA[0]) and olap > 0.75*(orfB[1]-orfB[0]):
            return True
    return False

    if orfA[0] == orfB[0] and orfA[1] == orfB[1]:
        return True
    elif orfA[0] == orfB[0]:
        if abs(orfA[1]-orfB[1]) < 0.1*min(orfA[1]-orfA[0], orfB[1]-orfB[0]):
            return True
    elif orfA[1] == orfB[1]:
        if abs(orfA[0]-orfB[0]) < 0.1*min(orfA[1]-orfA[0], orfB[1]-orfB[0]):
            return True
    return False

def compare(trueORFs, predORFs):
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

    '''
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
    #fn = len(true) - tp

    print '  Approved only:'
    print '    %% found: %0.2f (%d)' % (float(100*tp)/len(true), tp)
    print '    Genes missed: %d' % fn
    print '    Additional genes: %d' % fp
    '''

    
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

    print '  All predicted:'
    print '    %% found: %0.2f (%d)' % (float(100*tp)/len(true), tp)
    print '    Genes missed: %d' % fn
    print '    Additional genes: %d' % fp
    

