#!/usr/bin/env python
"""
dnolivieri: (updated: 18 jan 2016)
   - this is to identify the contigs that have some similarity
     to the 

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
from Bio.Blast.Applications import NcbitblastnCommandline as Tblastn
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

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

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

class CandidateContigs:
    def __init__(self, S):
        self.S = S
        self.in_query="./vrefs/consensus/consenso-V.fasta"


    def get_contigs_fromblastout(self, infile):

        def getList_fromblastOut(fbar):
            fp = open(fbar,"r")
            p=[]
            for lk in fp:
                s=lk.split()
                p.append( s[1].split("|")[3] )
            fp.close()
            return p

        sdir = os.path.dirname(infile)
        blast_List= list(find_files(sdir, "*blast_*") )
        #print blast_List
        p=[]
        for fb in blast_List: 
            print fb
            p.append(getList_fromblastOut(fb))


        #print "p=",p
        pbar = set(list(itertools.chain(*p)))
        p = list(pbar)
        p.sort()
        print "len(p)=",len(p)
        outfile = infile.replace(".fasta", "_cntg.txt")
        ofile = open(outfile, "w")
        for q in p: 
            ofile.write(q+"\n")

        ofile.close()




    def do_blast(self, in_file):
        tmpOutFile= in_file.replace(".fasta_", ".blast_")
        blastp_cmdline = Tblastn(query   = self.in_query,
                                subject = in_file,
                                evalue  = 1.e-3,   
                                outfmt=6,       
                                max_target_seqs=5000,
                                out=tmpOutFile)
        
        print blastp_cmdline
        stdout, stderr = blastp_cmdline()


    def preprocess_fasta(self,infile, totcnts):
        knt=0
        blck=0
        nncount =0
        rekcnt=0
        outfile = infile.replace(".fasta", ".fasta_"+str(blck) )            
        ofile = open(outfile, "w")
        for record in SeqIO.parse(infile, "fasta"):
            rekcnt+=1
            nncount+=len(record.seq)
            print knt, nncount
            if ((nncount > int(100.e6)) or (rekcnt>int(1.e5) )) and knt>0:
                ofile.close()
                print blck, knt
                print "Running tblast for it=",blck
                self.do_blast(outfile)

                blck+=1
                outfile = infile.replace(".fasta", ".fasta_"+str(blck) )            
                ofile = open(outfile, "w")
                nncount=0
                rekcnt=0

            SeqIO.write(record ,ofile, "fasta")
            knt+=1
        ofile.close()
        ## must do the last file
        self.do_blast(outfile)




    def run(self):
        pass
        """
        p = self.getList_fromblastOut(tmpOutFile) 
        pList = list(set(p))
        self.get_contigs_fromblastout(pList, in_file)
        #print "pList=", pList
        return pList
        """
            
## ---------------MAIN ----------------------------------
if __name__ == '__main__':
    parse_file = 'wgs_reptiles.json'
    json_data=open( parse_file )
    S = json.load(json_data)
    json_data.close()
    
    """
    mlist=["Ophiophagus_hannah","Ophisaurus_gracilis", 
           "Pogona_vitticeps_CEMB01"]
    """
    mlist=["Ophiophagus_hannah"]
    mlist=["Pogona_vitticeps_CEMB01"]
    C = CandidateContigs(S)
    for k in mlist:
        print k, S[k]["WGS"]
        #C.preprocess_fasta( S[k]["WGS"]   , 296399) 
        C.get_contigs_fromblastout(  S[k]["WGS"] )

        #C.do_blast( S[k]["WGS"] )
        

