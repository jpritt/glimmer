#! /usr/bin/env python
import sys
import himm
import himmTrain

def runHIMM(maxLength):
    trainer = himmTrain.HIMMTrain(maxLength)

    trainer.train()


    model = himm.HIMM(trainer.counts, trainer.singleProbs, maxLength)
    prediction = model.viterbi()

    with open('predictedMask.txt', 'w') as f:
        f.write(str(prediction) + '\n')