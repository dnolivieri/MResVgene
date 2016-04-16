#!/usr/bin/env python
"""
dnolivieri: (updated: 15 jan 2016)
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


class SelectRecords:
    def __init__(self, S):
        self.S = S

    def parse_fasta(self, inFile, locus):
        outFile = inFile.replace("_RF.fasta", "_RF_"+locus+".fasta")
        ofile = open( outFile,"w")        
        for record in SeqIO.parse(inFile, "fasta"):
            print record.id,
            if locus in record.id.split("|")[1]:
                print "----found"
                SeqIO.write(record ,ofile, "fasta")
            else:
                print 

        ofile.close()









    def run(self,  study, method, loci):
        self.parse_fasta ( self.S[study][method], loci )
        

            
## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    parse_file = 'parse_file.json'
    json_data=open( parse_file )
    S = json.load(json_data)
    json_data.close()
    
    method="rf"
    study="Pogona"
    loci="igk"  #  later this should be a list

    R = SelectRecords(S)
    R.run(study, method, loci)
