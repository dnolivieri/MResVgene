#!/usr/bin/env python
"""
   dnolivieri:  updated ...15 marzo 2015
      - get the regions from huge scaffolds.

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
from Bio.Alphabet import IUPAC

from scipy import *
import struct
import re
import json
import cPickle as pickle
from copy import deepcopy

import timeit
import operator



class getCntgRegions:
    def __init__(self, S):
        self.S = S

    def get_region(self, fbar): 
        outfile = fbar.replace(".fa", "_chr3.fasta")
        ofile = open(outfile, "w")

        Sbar = [ k for k in self.S.iterkeys() ]
        for record in SeqIO.parse(fbar, "fasta"):
            rec_desc = record.description
            cntg= rec_desc.split(":")[3]
            print cntg
            if cntg in Sbar: 
                print cntg, rec_desc
                for j in range(len(S[cntg]) ): 
                    q0 = S[cntg][j][0]
                    q1 = S[cntg][j][1]
                    print cntg, j, q0, q1
                    seq= record.seq[q0:q1]
                    recordB=SeqRecord(seq, id = record.id +"_"+str(j), description=record.description )
                    ofile.write(recordB.format("fasta"))

        ofile.close()

#-------------------
if __name__ == '__main__':

    infile="../VsGenomeB/mammals/Macaca_mulatta_ENSBL/Macaca_mulatta_MMUL1.fa"

    """
    S = { '1099214757507':[[0,5000]],
          '1099214148171':[[0,5000]],
          '1099548049584':[[0,1201553]],
          '1099214128018':[[0,30000]],
          '1099214732309':[[0,30000]],
          '10':[[65303593,66836272]],
          '13':[[88936957,90747293]],
          '13':[[112709235,113709678]],
          '7':[[83967455,85625679]],
          '3':[[179305555,180669600]]
          }
    """

    S = {
    '3':[[179305555,180669600]]
    }

    M = getCntgRegions(S)
    M.get_region(infile)


