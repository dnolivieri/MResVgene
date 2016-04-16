#!/usr/bin/env python
"""
   dnolivieri:  updated ...5-jan-2016
    -- re-write of the overlap algorithms.

"""

import collections
import numpy as np
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import struct
import re
from collections import defaultdict
from banyan import *
from copy import deepcopy
import timeit
import operator

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)



class Test:
    def __init__(self):
        pass


    def MRscore(self, pm):
        nlevel=2
        p = {}
        cnt=0
        for k in range(nlevel):
            p.update({k:[]})
            for m in range(np.power(2,k)):
                print k,m
                p[k].append(pm[cnt])
                cnt+=1
        print p

        sum_score = pm.sum()
        print "sum_score=", sum_score

        sbar=0.
        for k in range(1,nlevel):
            for m in range(np.power(2,k)):
                Delta = np.abs(p[0][0] - p[k][m])
                sigma=0.25
                sdelta = 1.0- np.exp( -np.power(Delta,2)/sigma ) 
                print k,m, p[k][m], Delta, sdelta
                sbar+= sdelta
                cnt+=1

        tot_score = sum_score - sbar 
        print "sum_score=", sum_score, "tot_score=", tot_score 



    def descriminate_by_exon(self, Q):
        #print "descriminate by exon:   Q=", Q
        ## Q[0] positions; Q[1] kmax; Q[2] max prob;  q[3] Pr_outclass
        s = [x[0] for x in Q]
        eprm = [x[1] for x in Q]
        score= [x[2] for x in Q]
        pclass = [ x[3] for x in Q] 
        sclass = [ x[4] for x in Q] 

        """
        print "-----------------"
        print "s=", s
        print "eprm=", eprm
        print "score=", score
        print "pclass=", pclass
        print "sclass=", sclass
        print "-----------------"
        """
        tree =  SortedSet( s, key=lambda (start, end): (end, (end - start)),updator=OverlappingIntervalsUpdator)
        r = set([])
        tbar=tree.copy()
        for sprm in s:
            a=tree.overlap(sprm)
            sindx=[s.index(x) for x in tree.overlap(sprm) ]
            z=[score[x] for x in sindx ]
            print 
            print "---sprm= ", sprm
            print "a=",a
            print "sindx=",sindx
            print "z=",z

            izmax = np.argmax(z)
            zmax  = np.max(z)
            print "izmax=", izmax
            print "zmax=", zmax
            
            for i in range(len(z)):
                if i!=izmax: 
                    tree.remove( s[sindx[i]] )
                    r.update([sindx[i]])


        r_list = [i for j, i in enumerate(Q) if j not in list(r)]
        print "before IN-CLASS"   
        print "r_list=", r_list
        print "rbar=",[x[0] for x in r_list]

        cand_list=r_list
        return cand_list







#-------------------
if __name__ == '__main__':


    Q= [((138384, 138697), 5, 1.8139954583665723, np.array([ 0.623,  0.535,  0.777]), 1.8139954583665723), 
        ((148946, 149248), 3, 2.0741487994214767, np.array([ 0.915,  0.819,  0.639]), 2.0741487994214767), 
        ((148955, 149248), 6, 1.5289169314516933, np.array([ 0.91 ,  0.559,  0.664]), 1.5289169314516933), 
        ((149507, 149802), 6, 2.0397966269792258, np.array([ 0.748,  0.674,  0.666]), 2.0397966269792258), 
        ((152502, 152790), 3, 1.8653537598968073, np.array([ 0.607,  0.763,  0.931]), 1.8653537598968073)]

    s=[ x[0] for x in Q]
    print s

    T = Test()
    cand_list = T.descriminate_by_exon(Q)
    for i in cand_list: 
        print i
    s=[ x[0] for x in cand_list]
    print s

    """

    T = Test()
    pm=np.array([ 0.854,  0.818,  0.76 ])
    T.MRscore(pm)
    """

