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
        outfile = fbar.replace(".fasta", "4.fasta")
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
                    recordB=SeqRecord(seq, id = cntg +"_"+str(j), description=record.description )
                    ofile.write(recordB.format("fasta"))

        ofile.close()

#-------------------
if __name__ == '__main__':


    ### For the reduced study of genes that didn't come out....17march


    infile="./MacacaMstudy/Macaca_mulatta_MMUL1_r.fasta"

    """
    S = { '1099214732309':[[1450,1850]],
          '7':[[1263900, 1264300], [987200, 987600], [376200, 376600], [839300, 839700] ],
          '3':[[849120,849520], [874900,875300]   ]
          }
    """

    S = {'7':[ [376200, 376600]]}


    M = getCntgRegions(S)
    M.get_region(infile)


    """
    def get_max( x ): 
        s = list(itertools.chain.from_iterable(x))
        s.sort()
        
        return np.array(s)

        
    res=[]


    S1=[[238,537]]
    res.append(  ('1099214757507', 0, 5000) )
    S2=[[797,1845]]
    res.append( ('1099214148171', 0, 5000) )
    S3=[[700393,701553], [654699,654992], [187218,187895] ]
    x = get_max(S3)
    res.append( ('1099548049584', 0, np.max(x)+500000) )
    print np.abs(np.min(x)-500000- np.max(x)+500000) 
    S4=[[627,1136]]
    res.append( ('1099214128018', 0, 30000) ) 

    S5=[[1261,1818]]
    res.append( ('1099214732309', 0,  30000) ) 

    Chr10=[ [66038886,66039290], [66222287,66222580], [66252992,66253527], [66335879,66336272], [66226801,66227112], 
            [65906638,65906946], [65929928,65930221], [65803593,65803889]]
    x = get_max(Chr10)
    res.append( ('10', np.min(x)-500000, np.max(x)+500000) )    
    print "chr10", np.abs(np.min(x)-500000- (np.max(x)+500000))


    Chr13=[[89436957,89437895],  [90246706,90247293], [89827964,89828257] ]
    x = get_max(Chr13)
    res.append( ('13', np.min(x)-500000, np.max(x)+500000) )    
    print "chr13", np.abs(np.min(x)-500000- (np.max(x)+500000) )


    Chr13B=[[113209235,113209678]]
    x = get_max(Chr13B)
    res.append( ('13B', np.min(x)-500000, np.max(x)+500000) )    
    print "chr13B", np.abs(np.min(x)-500000- (np.max(x)+500000) )


    Chr7=[ [84639642,84640144], [84695438,84695874],[85118821,85119249],[84991102,84991653],[84934218,84934823],
           [84496390,84496913],[84486442,84486950],[84467455,84468214],[84971778,84972336],[84716663,84717140],
           [84732052,84732598],[84564285,84564755],[85125178,85125679]]
    x = get_max(Chr7)
    res.append( ('7', np.min(x)-500000, np.max(x)+500000) )    
    print "chr7", np.abs(np.min(x)-500000- (np.max(x)+500000) )

  
    Chr3=[ [180086601,180112240],[179988975,179989609],[179951177,179951702],[179877783,179878247],
       [179805555,179806016],[179921910,179922345],[180168975,180169600],[180093458,180093897]]
    x = get_max(Chr3)
    res.append( ('3', np.min(x)-500000, np.max(x)+500000) )    
    print "chr3",  np.abs(np.min(x)-500000- (np.max(x)+500000) )


    print
    print res

    
    for k in res: 
        #print "'",k[0],"'",":[[",k[1], ",", k[2],"]]"
        print "'%s':[[%d,%d]]," % (k[0], k[1], k[2])


    """
