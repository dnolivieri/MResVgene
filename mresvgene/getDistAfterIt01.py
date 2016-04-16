#!/usr/bin/env python
"""
   dnolivieri: (11jan2016):
    - this does the iteration selection step for next iteration.
 """

import shutil
import numpy as np
import time
import bisect 
import os, fnmatch
import glob
import sys
import re

import itertools
import cPickle as pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Motif

from Bio import SeqFeature
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numpy.random import randn
import scipy

import KdeDist01 as KDE

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


class ProcessIteration: 
    def __init__(self, iterCnt, loci_classes):
        self.iterCnt = iterCnt
        self.inFile = "./bstrap/bstrp_iteration_"+str(self.iterCnt)+".fasta"
        self.loci_classes= loci_classes
        self.add_MotifSeqs()

    def add_MotifSeqs(self):
        text="inside add_Motif"
        self.m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
        self.m.add_instance(Seq("YYC",self.m.alphabet))
        self.m.add_instance(Seq("YFC",self.m.alphabet))
        self.m.add_instance(Seq("YLC",self.m.alphabet))
        self.m.add_instance(Seq("YIC",self.m.alphabet))
        self.m.add_instance(Seq("YHC",self.m.alphabet))
        self.m.add_instance(Seq("TFC",self.m.alphabet))

    def generate_training_set(self):
        pos_seqs = []

    def get_dist_data(self):
        B = {}
        D = {}
        for l in self.loci_classes:
            D.update({l:[]})
            B.update({l:0})


        for record in SeqIO.parse(self.inFile, "fasta"):
            rec_name=str(record.name)
            locus=rec_name.split("-")[2]
            rec_desc=str(record.description)
            mr_prob = float(rec_desc.split("|")[2])
            loci_prob = rec_desc.split("|")[3]
            probList=loci_prob.replace("(","").replace(")","").split(",")
            pbar = [ np.float(x) for x in probList ] 
            #print mr_prob, pbar
            p = record.seq

            if locus in self.loci_classes:
                D[locus].append( (mr_prob, pbar, rec_name) )

                z = p[33:45]
                z2 = p[18:30]
                z3 = p[-12:-1]
                if ('C' not in z2) and ('W' not in z):
                    B[locus]+=1

        return D, B





    def stats_sequences(self, D, X_Means):
        print "X_Means=X_means + Delta =", X_Means
        knt=0
        G = {}
        for l in self.loci_classes:
            G.update({l:[]})

        for k in self.loci_classes:
            print k, len(D[k]), X_Means[knt],
            jcnt=0
            for j in range(len(D[k])):
                #if (D[k][j][0] > X_Means[knt]) and (D[k][j][0] > 1.2):
                if (D[k][j][0] > X_Means[knt]):
                    #print j, D[k][j][2]
                    G[k].append(D[k][j][2])
                    jcnt+=1
            if len(D[k])>0:
                print "....jcnt=",jcnt,  float(jcnt)/float(len(D[k]))
            else:
                print "....jcnt=",jcnt,  float(jcnt)
            print
            knt+=1


    def select_sequences(self, D, X_Means):
        B={}
        for l in self.loci_classes:
            B.update({l:0})


        lfiles=[]
        lgfiles=[]
        for locus in self.loci_classes:
            f = "./bstrap/"+ locus + "_"+ str(self.iterCnt+1)+".fasta"
            g = "./bstrap/"+ locus + "_delta_"+ str(self.iterCnt+1)+".fasta"
            lfiles.append(open(f ,"w"))
            lgfiles.append(open(g ,"w"))

        outFile = self.inFile.replace(".fasta", "_r.fasta")
        ofile = open(outFile ,"w")
        knt=0
        for record in SeqIO.parse(self.inFile, "fasta"):
            rec_name=str(record.name)
            locus=rec_name.split("-")[2]
            rec_desc=str(record.description)
            mr_prob = float(rec_desc.split("|")[2])
            lindex = self.loci_classes.index( locus )

            #if (mr_prob > X_Means[lindex]) and (mr_prob > 1.2):
            p = record.seq
            z = p[33:45]
            z2 = p[18:30]
            z3 = p[-12:-1]

            cond1 =  ('C' in z2)
            cond2 =  ('W' in z) 
            cond3 =  len(list(self.m.search_instances(z3)))>0

            if (mr_prob > X_Means[lindex]) and cond1 and cond2 and cond3:
                SeqIO.write(record ,ofile, "fasta")
                SeqIO.write(record, lfiles[lindex], "fasta")
                SeqIO.write(record, lgfiles[lindex], "fasta")
    
                #and len(list(self.m.search_instances(z3)))>0:
                if ('C' not in z2) and ('W' not in z) and (('YYC' not in z3) 
                or ('YFC' not in z3) or ('YLC' not in z3) or ('YIC' not in z3) 
                or ('YHC' not in z3) or  ('TFC' not in z3)):

                    B[locus]+=1




                knt+=1
        print "Total=", knt
        print "B after=", B
        for k in self.loci_classes: 
            print k,  B[k]

        ofile.close()
        for k in range(len(lfiles)):
            lfiles[k].close()
            lgfiles[k].close()



    def run(self):
        D,B = self.get_dist_data()
        #print D


        make_plots=True
        Kd = KDE.KDEProbDistribution(D, self.loci_classes) 
        if make_plots:
            X_Means = Kd.get_kde_struct(show_plot=True)
            print "X_Means=",X_Means
        else:
            X_Means = Kd.get_kde_struct(show_plot=False)
            print "X_Means=",X_Means

            ## defines the training schedule
            """
            DeltaX = np.zeros(7)
            """
            #DeltaX = -0.1*np.ones(7) 
            DeltaX = 0.0*np.ones(7) 
            """
            DeltaX[0]=  -0.1  # ighv
            DeltaX[1]=  -0.1  # igkv
            DeltaX[2]=  -0.1 # iglv
            DeltaX[3]=  -0.1  # trav
            DeltaX[4]=  -0.1  # trbv
            DeltaX[5]=  -0.1  # trgv
            DeltaX[6]=  -0.1  # trdv
            """
            print "DeltaX=", DeltaX
            X_Means = X_Means + DeltaX


            X_Means = 2.0*np.ones(7)
            print "B=", B

            self.stats_sequences(D, X_Means)
            self.select_sequences(D, X_Means)


        
    def run_noKDE(self):
        D,B = self.get_dist_data()

        X_Means = 2.0*np.ones(7)

        X_Means[5] = 1.6
        X_Means[6] = 1.9

        #self.stats_sequences(D, X_Means)
        self.select_sequences(D, X_Means)





# ------------------------------------------
if __name__ == '__main__': 

    loci_classes = ["ighv", "igkv", "iglv", "trav", "trbv", "trgv", "trdv"]
    #loci_classes = ["ighv", "igkv", "iglv"]
    #loci_classes=["trav", "trbv", "trgv", "trdv"]
    n=4
    P= ProcessIteration(n,  loci_classes)
    P.run()
    #P.run_noKDE()
            
    """
    p, pos_threshold =  get_probs( str(record.description) )
    print "p=",p, "....cutoffs=",self.cutoffs,"...pos_threshold=", pos_threshold
    if sbar in good_seqs and pos_threshold:
    SeqIO.write(record ,self.ofile, "fasta")
    pos_seqs.append(sbar)
    """

    """
    exon1File = self.inFile.replace("outRF", "exon1")
    for record in SeqIO.parse(exon1File, "fasta"):
    sbar=str(record.name).split("exon1-")[1]
    if sbar in pos_seqs: 
    SeqIO.write(record ,self.nExon1, "fasta")
    """

    """
    Xcutoffs = Kd.get_kde_info()
    
    """
