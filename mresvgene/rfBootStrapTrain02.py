#!/usr/bin/env python
"""
  dnolivieri: (s)
   -- bootstrap training 

"""
import pylab as pl
import numpy as np
import sys

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cross_validation import cross_val_score
from numpy import genfromtxt, savetxt
import time
import itertools
import cPickle as pickle
import timeit
import pybedtools as pb

from Bio import SeqIO
from Bio.Seq import Seq

from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence

import json


AA = {1:'A', 2:'R',3:'N',4:'D',5:'C', 6:'Q', 7:'E', 8:'G', 9:'H', 10:'I',11:'L',  12:'K', 13:'M', 14:'F', 15:'P', 16:'S', 17:'T', 18:'W', 19:'Y',  20:'V',  21:'B', 22:'Z', 23:'X'}

rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

n_classes = 2
n_estimators = 1000
RANDOM_SEED = 13


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

# ------------------
class TrainEachExon:
    def __init__(self, S, Nestimators):
        self.S = S

        self.Nestimators = Nestimators
        self.make_training_matrices()


    def get_existing_signals(self, infile):
        data = np.load( infile )
        z = data['dp']
        return z


    def make_training_matrices(self):
        # background
        zB = self.get_existing_signals(self.S['Bkgnd']['bkg'] )
        print "zB.shape=", zB.shape
        valB = np.zeros( [zB.shape[0],1] )

        for exon in self.S['Exon'].iterkeys():
             zP = self.get_existing_signals(self.S['Exon'][exon] )
             valP = np.ones( [zP.shape[0],1] )
             print exon,  "zP.shape=", zP.shape
             train_signals = np.vstack([zP, zB])
             train_vals = np.vstack([valP, valB])
             outfile = "trainMat_"+ exon + ".pkl"

             A=TrainClassifier( train_signals, train_vals, self.Nestimators, outfile)


# ------------------
class TrainClassifier:
    def __init__(self, train_signals, train_vals, Nestimators, outfile, action=None):
        self.train_signals = train_signals
        self.Nestimators  = Nestimators
        self.outfile = outfile
        self.train_vals=  np.ravel(train_vals)

        print "do training"
        start_time = timeit.default_timer()
        rf=self.do_training()
        save_object(rf, self.outfile)
        elapsed = timeit.default_timer() - start_time
        print "ELAPSED=", elapsed

    def do_training(self):
        rf = RandomForestClassifier(n_estimators=self.Nestimators, oob_score=True)
        rf.fit(self.train_signals, self.train_vals)
        return rf

## ---------------MAIN ----------------------------------
if __name__ == '__main__':


    MHC1_Exon = 'MHC1_exon_npz.json'

    json_data=open( MHC1_Exon )
    S = json.load(json_data)
    json_data.close()

    T= TrainEachExon(S,  Nestimators=200)



