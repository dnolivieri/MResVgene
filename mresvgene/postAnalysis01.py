#!/usr/bin/env python
"""
   dnolivieri:  updated ...12 nov 2015
     - postAnalysis01.

     - Once the TransDecode/RSEM done, this 
       makes a summary.
      

      **  step 1:  obtain the 

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
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from scipy import *
import struct
import re
import json
import cPickle as pickle
from collections import defaultdict
import multiprocessing
from copy import deepcopy

import timeit
import operator

import errno

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


class PostAnalysis:
    def __init__(self, S):
        self.S = S
        self.consEvalue = 1.e-30
        #self.in_query="./data/IgM1_Pvitticeps.fasta"
        #self.in_query="./data/Vs_human.fasta"
        self.in_query="./gene_ref/homo_sapiens_IMGT.fasta"


    def getList_fromblastOut(self, infile):
        fp = open(infile,"r")
        p=[]
        for lk in fp:
            s=lk.split()
            p.append( (s[0],s[1]) )
        fp.close()
        return p



    def get_contigs_fromblastout(self, p, infile):
        print p
        igOut_file = infile.replace(".fasta.transdecoder.pep", "_transdecoder_Vs.fasta")
        pList = [ x[1] for x in p ]        
        ofile = open(igOut_file, "w")
        rec_cnt=0
        for record in SeqIO.parse(infile, "fasta"):
            if record.name in pList:
                SeqIO.write(record ,ofile, "fasta")
                rec_cnt+=1

        for record in SeqIO.parse(self.in_query, "fasta"):
            SeqIO.write(record ,ofile, "fasta")
            rec_cnt+=1                
        ofile.close()


        
        
    def do_blast(self, in_file):
        tmpOutFile= in_file.replace(".pep",".blast_out2")
        #in_query="./data/IgM1_Pvitticeps.fasta"

        blastp_cmdline = Blastp(query   = self.in_query,
                                subject = in_file,
                                evalue  = 1.e-20,   
                                #outfmt=2,       
                                max_target_seqs=10000,
                                out=tmpOutFile)
        
        #print blastp_cmdline
        stdout, stderr = blastp_cmdline()


        tmpOutFile= in_file.replace(".pep",".blast_out")
        blastp_cmdline = Blastp(query   = self.in_query,
                                subject = in_file,
                                evalue  = 1.e-20,   
                                outfmt=6,       
                                max_target_seqs=10000,
                                out=tmpOutFile)
        
        print blastp_cmdline
        stdout, stderr = blastp_cmdline()
        
        print "....Blastp output ......"
        fp = open(tmpOutFile, "r")
        for lp in fp: 
            print lp,


        p = self.getList_fromblastOut(tmpOutFile) 
        pList = list(set(p))
        self.get_contigs_fromblastout(pList, in_file)
        #print "pList=", pList
        return pList



    def getList_fromRSEM(self, infile, pList):
        fp = open(infile,"r")
        fpkm=[]
        pLabels=[ x[1].split("|")[0] for x in pList]
        #print "pLabels=", pLabels
        print "....RPKM from RSEM output ......"        
        for lk in fp:
            s=lk.split()[1].split(",")
            for x in s:
                if x in pLabels:
                    pindx= pLabels.index(x)
                    print pList[pindx][0], "...", x, ".....", lk.split()[6]
                    #print lk
                    fpkm.append( (pList[pindx][0], lk.split()[6] ) )

        fp.close()
        return fpkm
        


    def run(self):

        """
        sbar= S["diabetes1"]["t1"]
        sbar= S["preMs1"]["t1"]
        sbar= S["control1"]["t1"]
        print sbar
        p = self.do_blast(sbar["pep"])
        #print p
        r=self.getList_fromRSEM(sbar["rsem"], p)

        pklfile = sbar["pep"].replace(".pep", "_Vs.pkl")
        
        save_object(r, pklfile)
        """



        print "-------------------"
        for k in self.S.iterkeys():
            print "****", k
            for l in list(S[k]):
                print k, "...", l
                p = self.do_blast(S[k][l]["pep"])
                r=self.getList_fromRSEM(S[k][l]["rsem"], p)
                pklfile = S[k][l]["pep"].replace(".pep", "_Vs.pkl")
                save_object(r, pklfile)
                print "-------------------------------------------------"


        
# --------------------------------------
if __name__ == '__main__':

    RSEM_results = 'RSEM_summary.json'
    json_data=open(RSEM_results)
    S = json.load(json_data)
    json_data.close()


    R = PostAnalysis( S )
    R.run()



