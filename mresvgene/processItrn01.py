#!/usr/bin/env python
"""
   dnolivieri: (26dec2015):
    - this does the iteration selection step for next iteration.
    - adapted from cleanIteration
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
from Bio import SeqFeature
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numpy.random import randn
import scipy


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename



class KDEProbDistribution: 
    def __init__(self, cntIter):
        self.infile="./bstrap/bstrp_iteration_"+str(cntIter)+".fasta"
        self.cntIter = cntIter
        self.p=self.get_prob_fasta()



class ProcessIteration: 
    def __init__(self, p_iterN, loci_classes):
        self.p_iterN = p_iterN
        self.n_iterN = p_iterN + 1 
        self.loci_classes = loci_classes


    def getList_fromblastOut(self, infile):
        fp = open(infile,"r")
        p=[]
        for lk in fp:
            s=lk.split()
            if float(s[2])==100.:
                p.append( s[0] )
        fp.close()
        return p



    def get_blast_same(self, infile1, infile2):

        tmpOutFile= "blast_out"
        blastp_cmdline = Blastp(query   = infile2,
                                subject = infile1,
                                evalue  = 1.e-50,   
                                outfmt=6,       
                                max_target_seqs=100000,
                                out=tmpOutFile)
        
        print blastp_cmdline
        stdout, stderr = blastp_cmdline()
        
        print "....Blastp output ......"
        fp = open(tmpOutFile, "r")
        """
        for lp in fp: 
            print lp,
        """

        p = self.getList_fromblastOut(tmpOutFile) 
        pList = list(set(p))
        #self.get_contigs_fromblastout(pList, in_file)
        #print "pList=", pList
        return pList


    def only_difference_seqs(self):
        #2. Obtain the difference between previous iteration:
        if (self.p_iterN > 0):
            infile1="./bstrap/bstrp_iteration_"+str(self.p_iterN - 1 )+".fasta"
            infile2="./bstrap/bstrp_iteration_"+str(self.p_iterN)+".fasta"


            b1=[]
            for record in SeqIO.parse(infile1, "fasta"):     
                b1.append( record.id  )

            p1 = self.get_blast_same( infile1, infile2 )
            print "p1=",p1
            print "len(p1)", len(p1)

            bs1 = set(b1)
            ps1 = set(p1)
            s = bs1.difference(ps1)
            print s
            

    def get_new_sequences(self):
        fileList = list(find_files("./bstrap/", "*_"+str(self.p_iterN) +"_outRF.fasta"))
        #1. Get all sequences from this iteration and save in one file:
        outFile = "./bstrap/bstrp_iteration_"+str(self.p_iterN)+".fasta"
        ofile = open(outFile, "w")
        posCount=0
        for vFile in fileList:
            for record in SeqIO.parse(vFile, "fasta"):
                SeqIO.write(record, ofile, "fasta")
                posCount+=1

        ofile.close()
        print "posCount=", posCount

        #2 write each loci.
        for vFile in fileList:
            for loci in self.loci_classes: 
                outFile = vFile.replace(".fasta","_"+loci+".fasta")
                ofile = open(outFile, "w")
                for record in SeqIO.parse(vFile, "fasta"):
                    if record.id.split("-")[2]==loci: 
                        SeqIO.write(record, ofile, "fasta")
                ofile.close()


        # 3. first write the reference file: 
        for loci in self.loci_classes: 
            refFile = "./bstrap/"+loci + "_init.fasta"
            outFile = "./bstrap/"+loci +"_"+ str(self.p_iterN + 1) +".fasta"
            ofile = open(outFile, "w")
            for record in SeqIO.parse(refFile, "fasta"):
                SeqIO.write(record, ofile, "fasta")

            for vFile in fileList: 
                inFile = vFile.replace(".fasta","_"+loci+".fasta")                
                print loci, "inFile=", inFile 
                for record in SeqIO.parse(inFile, "fasta"):
                    SeqIO.write(record, ofile, "fasta")
            ofile.close()



# ------------------------------------------
if __name__ == '__main__': 

    presentN=1
    loci_classes=['ighv','igkv']
    P= ProcessIteration(presentN, loci_classes)
    P.get_new_sequences()
