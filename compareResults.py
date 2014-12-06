#! /usr/bin/env python

import sys
def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

#genome = readFASTA('~/Genomics/genomes/Haemophilus_influenzae/PittGG/sequence.fasta')
#genome = readFASTA(sys.argv[1])
#print genome[2915:4565]
#exit()

tp = 0
fp = 0
fn = 0

true = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        orf = line.rstrip().split('\t')
        true.append((int(orf[0]), int(orf[1])))

pred = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        orf = line.rstrip().split('\t')
        pred.append((int(orf[0]), int(orf[1])))

for orf in pred:
    found = False
    for o in true:
        if orf[1] == o[1]:# and orf[0] < o[0]:
            tp += 1
            found = True
            break
    if not found:
        fp += 1
for orf in true:
    found = False
    for o in pred:
        if orf[1] == o[1] and orf[0] < o[0]:
            found = True
            break
fn = len(true) - tp

print 'TP: %d' % tp
print 'FP: %d' % fp
print 'FN: %d' % fn