#!/usr/bin/env python
"""
   dnolivieri:  23 dec 2015

The Bootstrap supervised learning for VS  (adopted from the MHC work)

nohup ./bootstrap03.py glires 0 > bstrap/run_glires0_3.out &
nohup ./bootstrap03.py glires 1 > bstrap/run_glires1_3.out &
nohup ./bootstrap03.py primates 0 > bstrap/run_primates0_3.out &
nohup ./bootstrap03.py primates 1 > bstrap/run_primates1_3.out &

"""

import collections
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

import getFeatVec03 as FVec
import rfBootStrapTrain02 as RFtrain
import rfMapReduceMHC09 as RFpredict

class BootStrapML:
    def __init__(self, S, P, mlist ):
        self.S = S
        self.P = P
        self.mammalList = mlist
        self.Nestimators=500


    def get_sequences(self, fileList, iterCount):
        B =  FVec.GetFeatVectors(self.S,  desc_method='PDT2')
        posCount=0
        for exonFile in fileList:
            for record in SeqIO.parse(exonFile, "fasta"):
                posCount+=1

        for exonFile in fileList:
            D = B.descriptors_from_fasta(exonFile)
            npzoutfile = exonFile.replace("fasta","npz")
            np.savez( npzoutfile, dp=np.array(D) )

        
        bkgCount = posCount * 7
        print "posCount=", posCount,  "....bkgCount=", bkgCount
        zB = B.exon1_random_bkg(bkgCount)
        npzoutfile="./bstrap/bckgnd_exon1_"+str(iterCount)+ ".npz"
        np.savez(npzoutfile, dp=np.array(zB) )

        zB = B.exon2_random_bkg(bkgCount)
        npzoutfile="./bstrap/bckgnd_exon2_"+str(iterCount)+ ".npz"
        np.savez(npzoutfile, dp=np.array(zB) )

        zB = B.exon3_random_bkg(bkgCount)
        npzoutfile="./bstrap/bckgnd_exon3_"+str(iterCount)+ ".npz"
        np.savez(npzoutfile, dp=np.array(zB) )



    def make_training_matrices(self, fileList, iterCount):
        def get_signals(infile):
            data = np.load( infile )
            z = data['dp']
            return z


        exonList= ['exon1', 'exon2', 'exon3']
        for exon in exonList: 
            pos_file = "./bstrap/" +exon+ "_"+str(iterCount)+".npz"
            zP=get_signals(pos_file)
            valP = np.ones( [zP.shape[0],1] )

            bkg_file="./bstrap/bckgnd_"+exon+"_"+str(iterCount)+ ".npz"
            zB = get_signals(bkg_file)
            valB = np.zeros( [zB.shape[0],1] )

            print exon,  "zP.shape=", zP.shape, "zB.shape=", zB.shape
            train_signals = np.vstack([zP, zB])
            train_vals = np.vstack([valP, valB])

            outfile = "./bstrap/trainMat_"+ exon + "_"+str(iterCount)+ ".pkl"
            A= RFtrain.TrainClassifier( train_signals, train_vals, self.Nestimators, outfile)


    def TrainSubset(self, iterCount):
        #T = TrainEachExon(S,  Nestimators=200)
        pass

    def Predict(self, iterCount):

        mbar = []
        for m in self.mammalList: 
            mbar.append(m)
            M =  RFpredict.RunMHCpredict(self.P,  desc_method='PDT2',  mammalList= mbar )
            M.analyze_files(iterCount, 0.65)
            mbar = []


    def get_ResidualSet(self):
        
        def seqlist(fbar):
            seqList=[]
            for record in SeqIO.parse(fbar, "fasta"):
                seq = record.seq 
                seqList.append( (record.id, str(seq[0:12]) ) )
            return seqList

        z1 = seqlist(self.f1)
        z2 = seqlist(self.f2)
        y1 = [ x[1] for x in z1]
        y2 = [ x[1] for x in z2]

        r=[]
        for x in y2: 
            if x in y1: 
                i= y1.index(x)
                print i, z1[i]
                r.append(i)

        cand_list = [i for j, i in enumerate(z1) if j not in list(r)] 




    def run(self,n):
        flist=[]

        """
        for n in range(1):
            if n==0: 
                flist= [self.S[x] for x in self.S.iterkeys()  ]

        """
        flist= [self.S[x] for x in self.S.iterkeys() ]

        fbar = [ x.replace("init", str(n) )  for x in flist ]
        print n,  fbar
        
        #1. Sequences to Feature vecs
        self.get_sequences( fbar, n)
        
        #2. SL Training mats
        self.make_training_matrices( fbar, n)

        #3. Predict at iteration n
        self.Predict(n)
            
        #4. Determine new sequences
        





if __name__ == '__main__':
    import operator


    Vs_bstrap = 'mhc1_init_bstrap.json'
    json_data=open( Vs_bstrap )
    S = json.load(json_data)
    json_data.close()

    #MHC1_prediction = 'mhc1_mouse.json'
    #MHC1_prediction = 'mhc1_primates.json'


    family=sys.argv[1]
    hostmach = int(sys.argv[2])

    if family=='primates':
        MHC1_prediction = '../MHCwgscntgs/mhc_primates.json'

        if hostmach==0:
            mlist=["Aotus_JYKP01", "Callithrix_ACFV01", "Chlorocebus_AQIB01","Macacaf_AQIA01",
                   "Macacam_JSUE01","Macacan_JZLF01","Mandrillus_JYKQ01","Microcebus_ABDC02",
                   "Nomascus_ADFV01"]
        elif hostmach==1: 
            mlist=["Otolemur_AAQR03","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01",
               "Pongo_ABGA01","Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01",
               "Tarsius_ABRT02"]
        elif hostmach==2: 
            mlist=["TEST"]


    elif family=='glires':
        MHC1_prediction = '../MHCwgscntgs/mhc_glires.json'
        if hostmach==0:
            mlist=["Cavia_AAKN02", "Chinchilla_AGCD01", "Cricetulus_AMDS01", "Dipodomys_ABRO01", 
                   "Fukomys_AYUG01","Heterocephalus_AHKG01", "Jaculus_AKZC01"]
        elif hostmach==1: 
            mlist=["Mesocricetus_APMT01", "Microtus_AHZW01","Nannospalax_AXCS01","Octodon_AJSA01",
                   "Peromyscus_AYHN01", "Rattus_AABR07","Ictidomys_AGTP01"]
        elif hostmach==2: 
            mlist=["TEST"]



    else: 
        print "options are: [<primates> or <glires>]  <hostmach>"




    print mlist


    json_data2=open( MHC1_prediction )
    P = json.load(json_data2)
    json_data2.close()


    B = BootStrapML(S, P, mlist)
    n=0
    B.run(n)


    """
    for n in range(1):
        print n
        B.run(n)
    """


