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
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
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
        self.pdt3 = PDT3()
        self.Mstruct = MResStruct(self.nlevel)


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






    def hybrid_descriptors_from_sequence(self, seq):
        cnt=0
        D = {}
        for k in range(np.power(2,self.nlevel)-1):
            D.update( {k:[]} )

        descObject=PyPro.GetProDes(seq)
        seq_data = seq
        if ('X' not in seq_data ) and ('Z' not in seq_data) and ('B' not in seq_data):
            P  = self.Mstruct.get_pyramid(seq_data)
            knt=0
            for k in range(len(P)):
                qLambda = self.nlevel - k
                for kseq in P[k]:
                    T={}
                    if k==0: 
                        T = self.pdt3.get_freq_seq2vector(kseq)
                    else:
                        #T = self.pdt3.get_lambda_seq2vector( kseq, qLambda )
                        T = self.pdt3.get_freq_seq2vector(kseq)

                    Tx = [ T[str(x)]  for x in range(len(T)) ]
                    D[knt].append(Tx)
                    print cnt, k, knt, " ....", Tx[0:6],   kseq
                    knt+=1
                print
                cnt+=1
        return D




    def score_vgenes(self, speciesList): 
        self.rfmodels = self.get_models()

        for sbar in speciesList: 
            infile = S[sbar]["vref"]
            outfile=infile.replace("outV.fasta", "outV_mrfq.fasta")

            """
            infile = "./analysis/trees/all_S.fasta"
            outfile=infile.replace("S.fasta", "S_mr.fasta")
            """

            ofile = open(outfile, "w")

            for record in SeqIO.parse(infile, "fasta"):
                rec_name= record.name
                rec_seq=  record.seq.tostring()
                rec_desc= record.description

                Dbar = self.hybrid_descriptors_from_sequence(rec_seq)
                Y  =  np.array([Dbar])


                Ybar= []
                for k in range(np.power(2,self.nlevel)-1):
                    #Ybar.append(np.array( [list(itertools.chain(*x[k])) for x in Y ] ))
                    Ybar.append(np.array( [x[k] for x in Y ] ))
                
                ## could try to see if the vector length is >0

                LprobX=[]
                for loci in range(len(self.loci_classes)): 
                    probX=[]
                    for mrlevel in range(np.power(2,self.nlevel)-1):
                        try: 
                            probloci_level=self.rfmodels[loci][mrlevel].predict_proba(Ybar[:][mrlevel][0]) 
                            #print "loci=, mrlevel=, probloci_level=", loci, mrlevel, probloci_level
                            probX.append( probloci_level  )
                        except:
                            print "----ERROR in prob ----"
                            print rec_name, "   loci=",loci, " mrlevel=", mrlevel
                            print rec_seq
                            if len(Ybar[:][mrlevel])>0:
                                print "Ybar[0:10]=", Ybar[:][mrlevel][0]
                            else:
                                print "len(Ybar[:][mrlevel][0])= 0"
                            probloci_level = np.array( [[0., 0.]])
                            probX.append( probloci_level  )                            

                    LprobX.append(probX)


                #print "LprobX =", LprobX
                D={}
                for loci in self.loci_classes:
                    D.update({loci:[]})
                for lk in range(len(self.loci_classes)): 
                    loci = self.loci_classes[lk]
                    probX= LprobX[lk]
                    #print "probX=",probX

                    for j in range(probX[0].shape[0]):
                        Xbar = np.array([ probX[k][j][1] for k in range(len(probX)) ])
                        #print "Xbar=", Xbar
                        D[loci].append(Xbar)



                for j in range( len(D[self.loci_classes[0]])):  
                    q = [ D[x][j] for x in self.loci_classes ]
                    N_mrlevel = np.power(2,self.nlevel)-1
                    qMR = np.array([ self.MRscore(q[i])  for i  in range(len(self.loci_classes)) ])
                    iqMRmax = np.argmax(qMR)
                    qMRmax  = np.max(qMR)
                    qMx=q[iqMRmax]
                    print "j=", j,"  ", rec_name, "   ",  rec_seq, "   rho=", q, " iqMRmax=", iqMRmax, "   best_mrscores=",q[iqMRmax]
                    print "****** qMRmax=", qMRmax
                    print " ---------------------------------------"
                    
                    recordB=SeqRecord(record.seq, id = str(rec_name)+"-%1.3f[%1.3f-%1.3f-%1.3f]" % (qMRmax,qMx[0],qMx[1],qMx[2]), description="%1.3f"% qMRmax+"|"+str(iqMRmax) )
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
    mlist=["Chlorocebus_AQIB01","Gorilla_CABD02","Macacaf_AQIA01","Mandrillus_JYKQ01","Microcebus_ABDC01",
               "Nomascus_ADFV01","Panp_AJFE01","Pant_AACZ03","Papio_AHZZ01","Pongo_ABGA01", 
               "Propithecus_JZKE01","Rhinopithecus_JABR01","Saimiri_AGCE01","Tarsius_ABRT02"]
    mlist=["Tarsius_ABRT01"]

    #mlist=["Chlorocebus_AQIB01"]
    #mlist=["Macacaf_AQIA01"]
    #mlist=["Microcebus_ABDC01"]
    mlist=["Macacam_MMUL01"]

    V = VgeneDBPredict(S, 4, loci_classes)
    V.score_vgenes( mlist )


    """
    for sbar in speciesList: 
        infile = S[sbar]["vref"].replace(".fasta", "_mr.fasta")
        V.getKdedistributions(infile)
    """
