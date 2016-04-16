#!/usr/bin/env python
"""
   dnolivieri:  updated ...17 feb 2016
       - specially designed for looking at the VgeneDB sequences.
       - convert to feature vectors;  and give a prediction score.
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

import json
import cPickle as pickle
from collections import defaultdict
from banyan import *
import multiprocessing
from copy import deepcopy

import timeit
import operator


import  getMRfeatVec01 as FVec
from MResStruct02 import MResStruct
from PDT3method import PDT3


import KdeDist01 as KDE

# Note this redefinition to handle overlaps
#import regex as re
import re




class VgeneDBPredict:
    def __init__(self, S,  iterCount, loci_classes):
        self.S = S
        self.nlevel = 2
        self.iterCount = iterCount
        self.loci_classes = loci_classes
        self.rfmodels = []


    def get_models(self):
        rfmodels = []
        for loci in self.loci_classes: 
            nl=[]
            for mrlevel in range(np.power(2,self.nlevel)-1):
                matfile = "./bstrap/trainMat_" + loci + "_r"+str(mrlevel)+"_n"+ str(self.iterCount)+".pkl"
                print matfile
                fp = open(matfile, 'rb')
                nl.append( pickle.load(fp) )

            rfmodels.append( nl )
        return rfmodels

    def MRscore(self, pm):
        p = {}
        cnt=0
        for k in range(self.nlevel):
            p.update({k:[]})
            for m in range(np.power(2,k)):
                p[k].append(pm[cnt])
                cnt+=1
        sum_score = pm.sum()
        sigma=0.25
        sbar=0.
        for k in range(1,self.nlevel):
            for m in range(np.power(2,k)):
                Delta = np.abs(p[0][0] - p[k][m])
                sdelta = 1.0- np.exp( -np.power(Delta,2)/sigma ) 
                sbar+= sdelta
                cnt+=1

        tot_score = sum_score - sbar 
        return tot_score

    def get_sequence_buffer(self, infile):
        seqbuffer=[]
        for record in SeqIO.parse(infile, "fasta"):
            seqbuffer.append((record.id, record.seq))

        return seqbuffer

    def score_vgenes(self, speciesList): 
        self.rfmodels = self.get_models()
        F =  FVec.getMRFeatureVectors(self.S)
        for sbar in speciesList: 
            infile = S[sbar]["vref"]
            outfile=infile.replace("outV.fasta", "outV_mr.fasta")
            ofile = open(outfile, "w")
            seqbuffer= self.get_sequence_buffer(infile)

            Dbar=F.hybrid_descriptors_from_fasta(infile)
            print infile, " len(Dbar)=", len(Dbar),  len(Dbar[0])

            Y  =  np.array([Dbar])

            print "-----Step 2:  Reform matrix, Produce a numpy array for each MR level"
            Ybar= []
            for k in range(np.power(2,self.nlevel)-1):
                #Ybar.append(np.array( [list(itertools.chain(*x[k])) for x in Y ] ))
                Ybar.append(np.array( [x[k] for x in Y ] ))
                

            print "-----Step 3: For each MultRes level matrix, run prediction with approp model"
            LprobX=[]
            for loci in range(len(self.loci_classes)): 
                print "-----------------------------------"
                print "loci=", loci

                probX=[]
                for mrlevel in range(np.power(2,self.nlevel)-1):
                    try: 
                        probloci_level=self.rfmodels[loci][mrlevel].predict_proba(Ybar[:][mrlevel][0]) 
                        #print "loci=, mrlevel=, probloci_level=", loci, mrlevel, probloci_level
                        probX.append( probloci_level  )
                    except:
                        print "----ERROR in prob ----"
                        print "loci=",loci, " mrlevel=", mrlevel
                        if len(Ybar[:][mrlevel])>0:
                            print "Ybar[0:10]=", Ybar[:][mrlevel][0]
                        else:
                            print "len(Ybar[:][mrlevel][0])= 0"

                LprobX.append(probX)


            print "-----Step 4: Obtain Dict that summarizes data from each loci and all candidates"
            D={}
            for loci in self.loci_classes:
                D.update({loci:[]})
            for lk in range(len(self.loci_classes)): 
                loci = self.loci_classes[lk]
                probX= LprobX[lk]
                #print "probX[0].shape=", probX[0].shape
                #print "probX=", probX
                for j in range(probX[0].shape[0]):
                    # since proba is symmetric, only keep positive val...
                    Xbar = np.array([ probX[k][j][1] for k in range(len(probX)) ])
                    #print "Xbar=", Xbar
                    D[loci].append(Xbar)


            print "----Step 5: The logic for selecting Maximums:"


            for j in range( len(D[self.loci_classes[0]])):  
                q = [ D[x][j] for x in self.loci_classes ]
                N_mrlevel = np.power(2,self.nlevel)-1
                qMR = np.array([ self.MRscore(q[i])  for i  in range(len(self.loci_classes)) ])
                iqMRmax = np.argmax(qMR)
                qMRmax  = np.max(qMR)
                print "j=", j,"  ", seqbuffer[j][0], "   rho=", q, " iqMRmax=", iqMRmax
                print "****** qMRmax=", qMRmax
                print seqbuffer[j][1]
                recordB=SeqRecord(seqbuffer[j][1], id = str(seqbuffer[j][0])+"-%1.3f" % qMRmax, description="%1.3f"% qMRmax+"|"+str(iqMRmax) )
                ofile.write(recordB.format("fasta"))

            ofile.close()


                        
    def get_dist_data(self, infile):
        D = {}
        for l in self.loci_classes:
            D.update({l:[]})
        for record in SeqIO.parse(infile, "fasta"):
            rec_name=str(record.name)
            locus=rec_name.split("-")[2]
            mr_prob = float(rec_name.split("-")[3])
            if locus in self.loci_classes:
                D[locus].append( (mr_prob, rec_name) )
        return D



    def getKdedistributions(self, infile):
        D = self.get_dist_data(infile)
        make_plots=True
        Kd = KDE.KDEProbDistribution(D, self.loci_classes) 
        if make_plots:
            X_Means = Kd.get_kde_struct(show_plot=True)
            print "X_Means=",X_Means





#----------------------------------
if __name__ == '__main__':

    Vs_Loci = 'Vrefs_primates.json'

    json_data=open( Vs_Loci )
    S = json.load(json_data)
    json_data.close()

    loci_classes=[ 'ighv', 'igkv', 'iglv', 'trav','trbv','trgv', 'trdv']




    #speciesList=["Macacaf_AQIA01"]
    mlist=["Chlorocebus_AQIB01","Gorilla_CABD02","Macacaf_AQIA01",
               "Mandrillus_JYKQ01","Microcebus_ABDC02",
               "Nomascus_ADFV01","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01","Pongo_ABGA01", 
               "Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01","Tarsius_ABRT02"]

    mlist=["Chlorocebus_AQIB01"]
    V = VgeneDBPredict(S, 3, loci_classes)
    V.score_vgenes( mlist )


    """
    for sbar in speciesList: 
        infile = S[sbar]["vref"].replace(".fasta", "_mr.fasta")
        V.getKdedistributions(infile)
    """
