#!/usr/bin/env python
"""
   dnolivieri:  updated ...26 december 2015
      - modified to use the boostrap code.

      This does the prediction of Vs based 
      upon a multi-resolution training.

      - also should implement (eventually) the 
       bootstrap learning.
"""
from collections import defaultdict
import collections
import numpy as np
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
import struct
from copy import deepcopy
import re



def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)



def partition_func(s):
    D = {}
    for x in s:
        D.update({x[0]:[]})

    for x in s:
        D[x[0]].append(x[1])

    cand_list = []
    for k in D.iterkeys():
        print k, D[k]
        cand_list.append((k,D[k][0]))
    
    cand_list.sort()
    return cand_list




#----------------------------------
if __name__ == '__main__':

    """
    has this form:

    s= [((2423, 2738), Seq('SISILRRANQSIRVYF...FHG', ExtendedIUPACProtein())), 
    ((12429, 12736), Seq('SSGDIVMTQTP...QLP', ExtendedIUPACProtein())), 
    """
    #----before alg----
    Xbar= [ ((2423, 2738),'a'), ((12429, 12736),'b'), ((12621, 12911),'c'), 
            ((30471, 30771),'d'), ((30507, 30828),'e'), ((30528, 30828),'f'), 
            ((30537, 30828),'g'), ((39695, 39995),'h'), 
            ((39707, 39995),'i'), ((39695, 39995),'j'),  ((39695, 40005),'k'), 
            ((39707, 39995),'l'), ((39707, 40005),'m'), ((42076, 42368),'n'), 
            ((51597, 51893),'p')]


    """
    sbar=[ x[0] for x in Xbar]
    for x in sbar: 
        if x[0] > 30000 and x[0] < 40000:
            print x

    """


    X = partition_func(Xbar)

    print X

    """
    print 
    print X

    for x in X: 
        if x[0] > 30000 and x[0] < 40000:
            print x


    sx=[(30471, 30771),(30507, 30828),(30528, 30828),(30537, 30828),(39695, 39995),(39707, 39995),(39695, 39995),(39695, 40005),(39707, 39995),(39707, 40005)]

    print sx
    sxprm=list(set(sx))
    sxprm.sort()
    print sxprm

    tx=[(30471, 30771),(30507, 30828),(30528, 30828),(30537, 30828),(39695, 40005),(39707, 40005)]
    print tx
    """
