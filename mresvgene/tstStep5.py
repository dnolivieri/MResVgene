#!/usr/bin/env python
"""
   dnolivieri:  updated ...8 december 2015
      vregMRmodel:
      The specific multi-resolution (MR) model for V-regions.
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


nlevel=2


def get_MRseqs(loci_classes, D):
    print "STEP 5 ---------------------------"
    seqs=[]
    for j in range( len(D[loci_classes[0]])):  
        rho = [ D[x][j] for x in D.iterkeys() ]
        tau = np.array([ D[x][j].sum() for x in D.iterkeys() ])
        loci_kmax = np.argmax(tau)
        loci_ymax = np.max(tau)
        if ( np.max(rho[loci_kmax])  > 0.4 ):
            print j, "....rho=", rho, "    tau=", tau, ".....", loci_kmax, loci_ymax,"***",rho[loci_kmax], np.max(rho[loci_kmax])
            rx=rho[loci_kmax]
            sbar=True
            for k in range(1,np.power(2,nlevel)-1):
                sbar = sbar and (np.abs(rx[k]- rx[0]) < 0.2)
                print "     ",k, rx[0], rx[k],"...", np.abs(rx[k]- rx[0]), sbar

            if sbar==True:
                seqs.append(j)
        else:
            print

    return seqs



    


if __name__ == '__main__':
    D= {'ighv': [array([ 0.7  ,  0.278,  0.49 ]), array([ 0.468,  0.284,  0.174]), 
                 array([ 0.31 ,  0.224,  0.016]), array([ 0.33 ,  0.198,  0.018]), 
                 array([ 0.312,  0.034,  0.094]), array([ 0.742,  0.632,  0.654]), 
                 array([ 0.82 ,  0.838,  0.916]), array([ 0.38 ,  0.606,  0.494]), 
                 array([ 0.882,  0.36 ,  0.216]), array([ 0.876,  0.564,  0.208]), 
                 array([ 1.   ,  0.998,  0.996]), array([ 0.034,  0.   ,  0.002]), 
                 array([ 0.026,  0.   ,  0.004]), array([ 0.026,  0.002,  0.002]), 
                 array([ 0.01,  0.  ,  0.  ]), array([ 0.078,  0.008,  0.058]), 
                 array([ 0.08 ,  0.006,  0.076]), array([ 0.076,  0.   ,  0.082]), 
                 array([ 0.064,  0.   ,  0.058]), array([ 0.05 ,  0.   ,  0.038])], 
        'igkv': [array([ 0.078,  0.06 ,  0.132]), array([ 0.182,  0.034,  0.018]), 
                 array([ 0.152,  0.028,  0.012]), array([ 0.12 ,  0.032,  0.006]), 
                 array([ 0.026,  0.058,  0.034]), array([ 0.096,  0.166,  0.05 ]), 
                 array([ 0.068,  0.088,  0.038]), array([ 0.092,  0.05 ,  0.058]), 
                 array([ 0.28 ,  0.068,  0.052]), array([ 0.266,  0.066,  0.058]), 
                 array([ 0.312,  0.14 ,  0.156]), array([ 0.018,  0.008,  0.004]), 
                 array([ 0.018,  0.016,  0.004]), array([ 0.024,  0.026,  0.004]), 
                 array([ 0.016,  0.018,  0.002]), array([ 0.01 ,  0.014,  0.01 ]), 
                 array([ 0.012,  0.014,  0.01 ]), array([ 0.008,  0.014,  0.008]), 
                 array([ 0.014,  0.   ,  0.012]), array([ 0.008,  0.   ,  0.006])]}


    loci_classes=['ighv', 'igkv']
    seqs = get_MRseqs(loci_classes, D)


    print "seqs=",seqs
