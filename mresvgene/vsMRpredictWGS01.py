#!/usr/bin/env python
"""
   dnolivieri:  23 dec 2015
       Instead of bootstrapping, this 
       just does the prediction using matrices

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

import  mrvPredict02WGS as RFpredict

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class vsMRpredictML:
    def __init__(self, loci_classes, P, mlist):
        self.P = P
        self.loci_classes = loci_classes
        self.speciesList = mlist
        self.Nestimators=250
        #self.bkg_total = 10000
        self.bkg_total = 6500


    def Predict(self, iterCount):
        mbar = []
        threshold = 0.5
        """
        cntgs=['chromosome:MMUL_1:3:1:196418989:1']
        cntgs=['scaffold:MMUL_1:1099548049584:264277:477440:1']
        """

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
        print "#3._____ Predict at iteration n"
        self.Predict(n)

        
# -----------------------------------------------
if __name__ == '__main__':


    Vs_wgs_prediction = 'wgs_primates.json'
    #mlist=["Macacam_MMUL1"]
    mlist=["Macacam_MMUL1S"]
    print mlist

    json_data2=open( Vs_wgs_prediction )
    P = json.load(json_data2)
    json_data2.close()
    loci_classes = ["ighv", "igkv", "iglv", "trav", "trbv", "trgv", "trdv"]


    B = vsMRpredictML(loci_classes, P, mlist)


    n=4
    print "**** N = ", n, "  ****"
    B.run(n)
