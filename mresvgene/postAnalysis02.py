#!/usr/bin/env python
"""
   dnolivieri:  updated ...12 nov 2015
     - postAnalysis01.

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


class KDEProbDistribution: 
    def __init__(self, cntIter):
        self.infile="./bstrap/bstrp_iteration_"+str(cntIter)+".fasta"
        self.cntIter = cntIter
        self.p=self.get_prob_fasta()


class PostAnalysis:
    def __init__(self, S):
        self.S = S
        self.consEvalue = 1.e-30


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




    def run(self):
        
# --------------------------------------
if __name__ == '__main__':

    RSEM_results = 'RSEM_summary.json'
    json_data=open(RSEM_results)
    S = json.load(json_data)
    json_data.close()


    R = PostAnalysis( S )
    R.run()



