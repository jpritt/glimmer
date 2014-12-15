#!/usr/bin/env python

import sys
from collections import defaultdict

def readFASTA(filename):
    seqs = dict()
    with open(filename, 'r') as f:
    	label = ''
    	temp_seq = ''
        for line in f:
            if line[0] == '>':
            	if len(label) > 0:
            		seqs[label] = temp_seq
            		temp_seq = ''
            	label = line[1:].rstrip()
            else:
            	temp_seq += line.rstrip()
    return seqs

def readExonDef(filename):
    seqs = defaultdict(list)
    with open(filename, 'r') as f:
    	label = ''
    	temp_seq = ''
        for line in f:
    		exon_loc = line.rstrip().split(' ')
    		if len(exon_loc) == 3:
	    		label = exon_loc[0]
	    		start = exon_loc[1]
	    		end = exon_loc[2]
	    		seqs[label].append((start, end))
    return seqs

seqs = readFASTA(sys.argv[1])

print seqs['X89714']

coding_seqs = readExonDef(sys.argv[2])

print coding_seqs['X89714']
