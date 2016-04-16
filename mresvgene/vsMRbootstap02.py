#!/usr/bin/env python
"""
   dnolivieri:  23 dec 2015
     Bootstrap workflow for the MR ensemble RF code.

"""

import collections
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json
import cPickle as pickle
from collections import defaultdict
from banyan import *
import multiprocessing
from copy import deepcopy

import timeit
import operator

import  getMRfeatVec01 as FVec
#import  mrvTrain01 as RFtrain
import  mrvTrainSG02 as RFtrain
import  mrvPredict03 as RFpredict
import  processItrn02  as pItrn

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class vsBootStrapML:
    def __init__(self, S, loci_classes, P, mlist):
        self.S = S
        self.P = P
        self.loci_classes = loci_classes
        self.speciesList = mlist
        self.Nestimators=250
        #self.bkg_total = 10000
        self.bkg_total = 6500

    def convert_full_background(self):
        F =  FVec.getMRFeatureVectors(self.S)
        #bkgndFile = "./dataVgeneDB/noMotifsBcknd.fasta"
        #bkgndFile = "./dataVgeneDB/bckgnd_run_r6500.fasta"
        bkgndFile = "./dataVgeneDB/All_bkg_m.fasta"
        #Zbar = F.background_from_fasta(bkgndFile, self.bkg_total)
        ## Try the hybrid method
        Zbar = F.hybrid_descriptors_from_fasta(bkgndFile, self.bkg_total)

        pkloutfile="./bstrap/bckgnd_n" + str(self.bkg_total) + ".pkl"
        save_object(Zbar, pkloutfile)

    def get_existing_data(self, infile):
        qp = open(infile, 'rb')
        D = pickle.load(qp)
        return D


    def convert_sequences(self, iterCount):
        print "-----inside convert_sequences------"
        F =  FVec.getMRFeatureVectors(self.S)
        G = {}
        for loci in self.loci_classes:
            posCount=0
            G.update({loci:0})

            """
            if iterCount==0:
                lociFile = "./bstrap/" + loci + "_" + str(iterCount) +".fasta"
            else:
                lociFile = "./bstrap/" + loci + "_delta_" + str(iterCount) +".fasta"
            """
            ## it will always be the first one!... I still have not done the delta correctly
            lociFile = "./bstrap/" + loci + "_0.fasta"

            for record in SeqIO.parse(lociFile, "fasta"):
                posCount+=1

            print "lociFile=", lociFile,  " posCount=", posCount

            if iterCount>0:
                prev_loci_file = "./bstrap/" + loci + "_" + str(iterCount - 1 ) +".pkl"
                p = self.get_existing_data(prev_loci_file)
                zp= np.array(p[0])
                posCount += zp.shape[0]
                print "zp.shape[0]=", zp.shape[0], "posCount=", posCount



            if posCount * 3 < 10000:
                G[loci] = posCount * 3
            else:
                G[loci] = 9999
            """
            if G[loci] < 200:
                G[loci] = 200
            """
                
        print "G=", G
        print "---getting sequences from background ---"
        bkgndFile = "./bstrap/bckgnd_n" + str(self.bkg_total) + ".pkl"
        zB = self.get_existing_data(bkgndFile)
        print "---starting loci conversion ---"
        for loci in self.loci_classes:
            lociFile = "./bstrap/" + loci + "_" + str(iterCount) +".fasta"
            #Dbar=F.descriptors_from_fasta(lociFile)
            ## try hybrid method
            Dbar=F.hybrid_descriptors_from_fasta(lociFile)
            print loci, " len(Dbar)=", len(Dbar),  len(Dbar[0])

            if iterCount>0:
                prev_loci_file = "./bstrap/" + loci + "_" + str(iterCount - 1 ) +".pkl"
                print "prev_loci_file=", prev_loci_file
                yP = self.get_existing_data(prev_loci_file)
                if len(Dbar)>0:
                    ybar={}
                    print "updating "
                    for k in range(len(yP)):
                        ybar.update({k:yP[k]+Dbar[k]})
                        print len(yP[k]), len(Dbar[k]), len(ybar[k])
                    yP = ybar

            else:
                yP = Dbar

            for k in range(len(yP)):
                print len(yP[k])
            
            pkloutfile=lociFile.replace(".fasta", ".pkl")
            print "pkloutfile=", pkloutfile
            save_object(yP, pkloutfile)

            rand_samp=True
            zbar = {}
            jz=rnd.sample(range(0,self.bkg_total),int(G[loci]))
            for k in range(len(zB)):
                # get random samples...
                print G[loci]
                if rand_samp:
                    #zbar.update({k:zB[k][jz] } )
                    zbar.update({k:[zB[k][zkk] for zkk in jz]} )
                else:
                    zbar.update({k:zB[k][0:G[loci]] } )

                x= np.array(zbar[k])
                print x.shape 

            pkloutfile="./bstrap/bckgnd_"+loci+"_"+ str(iterCount)+ ".pkl"
            save_object(zbar, pkloutfile)



    def make_training_matrices (self, iterCount ): 
        print "---making loci based training matrices ----"
        for loci in self.loci_classes:
            loci_file = "./bstrap/" + loci + "_" + str(iterCount) +".pkl"
            print "loci_file=", loci_file
            D = self.get_existing_data( loci_file )
            bkg_file="./bstrap/bckgnd_"+loci+"_"+ str(iterCount)+ ".pkl"
            B = self.get_existing_data( bkg_file )

            for k in range(len(D)):
                zP=np.array(D[k])                
                zB=np.array(B[k])
                valP = np.ones( [zP.shape[0],1] )
                valB = np.zeros( [zB.shape[0],1] )
                print k, zP.shape, zB.shape, valP.shape, valB.shape
                train_signals = np.vstack([zP, zB])
                train_vals = np.vstack([valP, valB])
                print "train_signals.shape=", train_signals.shape 
                print "train_vals.shape=", train_vals.shape
                out1="./bstrap/train_signal_"+ loci + "_"+ str(k) +".pkl"
                out2="./bstrap/train_val_"+ loci + "_"+ str(k) +".pkl"
                outfile = "./bstrap/trainMat_" + loci + "_r"+str(k) +"_n"+ str(iterCount)+".pkl"
                print outfile
                A= RFtrain.TrainClassifier( train_signals, train_vals,  outfile)




    def Predict(self, iterCount):
        mbar = []
        #threshold = 0.95 - iterCount*0.1
        ##  use the 0.5 threshold to get a distribution;  the iteration must be done by hand....
        threshold = 0.5
        #cntgs=['AQIA01064159.1']
        #cntgs=['AQIB01146122.1']
        """
        cntgs=['CABD02231421.1']
        cntgs=['CABD02231427.1']
        cntgs=['AJFE01012238.1']
        cntgs=['AHZZ01042080.1']
        cntgs=['ABGA01394058.1']
        cntgs=['ABGA01252219.1']
        cntgs=['ABGA01275597.1']
        cntgs=['JZKE01249934.1']
        cntgs=['JABR01071860.1']
        cntgs=['JABR01020049.1']
        cntgs=['JABR01016432.1']
        cntgs=['ABRT010282393.1']
        """
        #cntgs=['ABRT010078850.1']

        #cntgs=['chromosome:MMUL_1:3:1:196418989:1']
        cntgs=None
        for m in self.speciesList: 
            mbar.append(m)
            M =  RFpredict.RunVregPredict(self.P,  speciesList= mbar, check_contigs=True, contigs=cntgs)
            M.analyze_files(iterCount, self.loci_classes, threshold)
            mbar = []


    def PrepareNext(self, iterCount ):
        P= pItrn.ProcessIteration(iterCount, self.loci_classes)
        P.run()

    def run(self,n):
        # 0. Do the background once: 
        #self.convert_full_background()
        
        """
        print "#1.____ Sequences to Feature vecs"
        self.convert_sequences( n )

        print "#2._____ SL Training mats"
        self.make_training_matrices( n )
        """

        print "#3._____ Predict at iteration n"
        self.Predict(n)

        """
        print "#4.____ Determine new sequences"
        self.PrepareNext( n )
        """
        
# -----------------------------------------------
if __name__ == '__main__':

    Vs_bstrap = 'Vs_init_bstrap_ighv.json'
    json_data=open( Vs_bstrap )
    S = json.load(json_data)
    json_data.close()

    #family=sys.argv[1]
    #hostmach = int(sys.argv[2])
    family = 'primates'
    #family = 'reptiles'

    if family=='primates':
        #Vs_wgs_prediction = 'wgs_primates_ig.json'
        ### this is the full 
        #Vs_wgs_prediction = 'wgs_primates_tr.json'
        Vs_wgs_prediction = 'wgs_primates.json'


        mlist=["Chlorocebus_AQIB01","Gorilla_CABD0","Macacaf_AQIA01",
               "Macacam_AANU01","Macacan_JZLF01","Mandrillus_JYKQ01","Microcebus_ABDC02",
               "Nomascus_ADFV01","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01","Pongo_ABGA01", 
               "Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01","Tarsius_ABRT02"]

        """
        mlist=["Chlorocebus_AQIB01","Gorilla_CABD0","Macacaf_AQIA01",
               "Mandrillus_JYKQ01","Microcebus_ABDC02",
               "Nomascus_ADFV01","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01","Pongo_ABGA01", 
               "Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01","Tarsius_ABRT02"]
        """
        #mlist=["Chlorocebus_AQIB01"]
        #mlist=["Macacam_AANU01", "Mandrillus_JYKQ01"]
        #mlist=["Macacam_AANU01"]
        mlist=["Chlorocebus_AQIB01","Gorilla_CABD02","Macacaf_AQIA01",
               "Mandrillus_JYKQ01","Microcebus_ABDC01",
               "Nomascus_ADFV01","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01","Pongo_ABGA01", 
               "Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01","Tarsius_ABRT01"]

        #mlist=["Tarsius_ABRT01"]

        #mlist=["Macacam_MMUL1"]
    elif family=='glires':
        Vs_wgs_prediction = 'wgs_glires.json'
        mlist=[ ]


    elif family=='reptiles':
        Vs_wgs_prediction = 'wgs_reptiles.json'
        mlist=["Ophiophagus_hannah"]        


    print mlist

    json_data2=open( Vs_wgs_prediction )
    P = json.load(json_data2)
    json_data2.close()


    #loci_classes = ["ighv", "igkv", "iglv"]

    #loci_classes = ["trav", "trbv", "trgv", "trdv"]
    loci_classes = ["ighv", "igkv", "iglv", "trav", "trbv", "trgv", "trdv"]
    B = vsBootStrapML(S, loci_classes, P, mlist)


    """
    for n in range(0,3):
        print "**** N = ", n, "  ****"
        B.run(n)
    """



    n=4
    print "**** N = ", n, "  ****"
    B.run(n)
