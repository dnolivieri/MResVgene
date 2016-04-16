#!/usr/bin/env python
"""
   dnolivieri:  updated ...7 dec 2015
      - based on my previous code, but this does the pyramid.
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
from Bio import SeqIO
from Bio.Seq import Seq

from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json

import genPosMRVecs02 as gPosV

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class MResTraining:
    def __init__(self, S, loci_classes, pos_exists=False, neg_exists=False):
        self.S = S
        self.loci_classes  = loci_classes
        self.pos_exists= pos_exists
        self.neg_exists= neg_exists
        self.posFile    = 'train_posSignal'
        self.bckgndFile = 'train_bckgndSignal'
        self.plot_featurevecs = False
        self.form_MRtrainigvectors()

        """
        train_signals, train_vals =  self.form_featurevecs()
        #train_signals, train_vals =  self.form_featurevecs_fromfiles()
        npzsignals_out = "train_signals.npz"
        np.savez(npzsignals_out, dp=train_signals )
        npzvals_out = "train_vals.npz"
        np.savez(npzvals_out, dp=train_vals )

        do_training=True
        if do_training:
            A=TrainClassifier( train_signals, train_vals, 1000 )

        """
    def get_existing_data(self, infile):
        qp = open(infile, 'rb')
        D = pickle.load(qp)
        return D

    def form_MRtrainigvectors(self):
        # background
        bkg_file=self.S['All']['bkg1'].replace(".fasta",".pkl")
        print "background=", bkg_file 
        zB1 = self.get_existing_data( bkg_file )

        n_iter=0
        print "getting loci data"
        for loci in self.loci_classes:
            loci_file = self.S['All'][loci].replace(".fasta",".pkl")
            print "loci_file=", loci_file
            D = self.get_existing_data( loci_file )
            print "len(D)=", len(D)

            for k in range(len(D)):
                zP=np.array(D[k])                
                zB=np.array(zB1[k])
                valP = np.ones( [zP.shape[0],1] )
                valB = np.zeros( [zB.shape[0],1] )
                print k, zP.shape, zB.shape, valP.shape, valB.shape
                train_signals = np.vstack([zP, zB])
                train_vals = np.vstack([valP, valB])
                print "train_signals.shape=", train_signals.shape 
                print "train_vals.shape=", train_vals.shape
                out1="./bstrap/train_signal_"+ loci + "_"+ str(k) +".pkl"
                out2="./bstrap/train_val_"+ loci + "_"+ str(k) +".pkl"
                #save_object(train_signals, out1)
                #save_object(train_vals, out2)

                #outfile = "./bstrap/trainMat_"+ loci + "_"+ str(k) +".pkl"
                outfile = "./bstrap/trainMat_" + loci + "_r"+str(k) +"_n"+ str(n_iter)+".pkl"

                print outfile
                A=TrainClassifier( train_signals, train_vals,  outfile)




# ------------------
class TrainClassifier:
    def __init__(self, train_signals, train_vals,  outfile):
        self.train_signals = train_signals
        self.train_vals    = train_vals
        self.Nestimators  = 500
        self.outfile = outfile

        #self.train_vals=  train_vals.reshape( [train_vals.shape[0], 1])
        self.train_vals=  np.ravel(train_vals)

        do_training=True
        if do_training:
            print "do training"
            start_time = timeit.default_timer()
            rf=self.do_training()
            self.get_cross_val_score()
            #save_object(rf, r'train_Matrix.pkl')
            save_object(rf, self.outfile)
            elapsed = timeit.default_timer() - start_time
            print "ELAPSED=", elapsed

        count_vectors=True
        if count_vectors:
            self.obtain_trainig_set_size()

    def do_training(self):
        rf = RandomForestClassifier(n_estimators=self.Nestimators, oob_score=True)
        #rf = ExtraTreesClassifier(n_estimators=self.Nestimators, bootstrap=True, oob_score=True)
        rf.fit(self.train_signals, self.train_vals)

        return rf


    def obtain_trainig_set_size(self):
        print np.count_nonzero(self.train_vals)
        sbar=set(np.ravel(self.train_vals).astype(int))

        for k in list(sbar):
            print k, np.extract( self.train_vals==k, self.train_vals).size

    def get_feature_importance(self, rf):
        impList=[]
        importances = rf.feature_importances_
        std = np.std([tree.feature_importances_ for tree in rfmodel.estimators_], axis=0)
        indices = np.argsort(importances)[::-1]
        print "indices.shape=", indices.shape

        #Print the feature ranking
        print("Feature ranking:")            
        xbar=[]
        ybar=[]
        for f in range(5):
            print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
            #xbar.append(  indices[f] )
            #ybar.append(importances[indices[f]] )


    def get_cross_val_score(self):
        rf = RandomForestClassifier(n_estimators=self.Nestimators, oob_score=True)
        scores = cross_val_score(rf,  self.train_signals,  self.train_vals , cv=5)
        print "scores=", scores, 
        print "mean score=", scores.mean()



## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    Vs_Loci = 'Vs_Loci_fasta.json'

    json_data=open( Vs_Loci )
    S = json.load(json_data)
    json_data.close()

    #loci_classes=[ 'ighv', 'iglv', 'igkv', 'trav','trbv','trgv', 'trdv']
    loci_classes=[ 'ighv', 'igkv' ]


    MResTraining( S, loci_classes )
