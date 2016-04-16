#!/usr/bin/env python
"""
   dnolivieri:  updated ...8 december 2015
      vregMRmodel

      The specific multi-resolution (MR) model for V-regions.
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
            b=[s.index(x) for x in tree.overlap(sprm) ]
            c=[eprm[x] for x in b ]
            d=[score[x] for x in b ]

            pbar = [pclass[x] for x in b]
            spbar = [sclass[x] for x in b]

            gX=sorted(list_duplicates(c))
            if len(gX)>0:
                g=sorted(list_duplicates(c))[0][1]
                z=[d[i] for i in g ]
                z=   [d[i] for i in g ]
                zp = [pbar[i] for i in g]
                zs = [spbar[i] for i in g]
                izmax = self.maxLk_interval( z, zs )

                for i in g:
                    #if i!= g[np.argmax(z)]:
                    if i!= g[izmax]: 
                        tree.remove( s[b[i]] )
                        r.update([b[i]])
                else:
                    z=[]

        r_list = [i for j, i in enumerate(Q) if j not in list(r)]
        c=[ x[1] for x in r_list ]

        print "before IN-CLASS"   
        print "r_list=", r_list
        print "c=", c

        # ---- IN CLASS Overlaps ----------------

        tmpbuff=[]
        rindx_buffer=[]
        for i in range(1,len(c)):
            if (i-1) not in tmpbuff:
                tmpbuff.append(i-1)
            if c[i]== c[i-1]:
                if i not in tmpbuff:
                    tmpbuff.append(i)
            else:
                if len(tmpbuff)>0:
                    z = [ r_list[x][2] for x in tmpbuff ]
                    imax=np.argmax(z)
                    jmax=tmpbuff[imax]
                    for j in tmpbuff:
                        if j!=jmax:
                            rindx_buffer.append(j)

                tmpbuff=[]
            lastI=i

        merged = rindx_buffer
        Q = [i for j, i in enumerate(r_list) if j not in merged]
        print "after IN-CLASS   Q =", Q

        ## ---- Out-Of-Class Overlaps -------------------------

        s = [ x[0] for x in Q ]
        eprm = [x[1] for x in Q]
        score= [x[2] for x in Q]

        tree =  SortedSet( s, key=lambda (start, end): (end, (end - start)),updator=OverlappingIntervalsUpdator)
        r = set([])
        for sprm in s:
            a=tree.overlap(sprm)
            b=[s.index(x) for x in tree.overlap(sprm) ]
            c=[eprm[x] for x in b ]
            d=[score[x] for x in b ]

            if len(b)>1:
                for i in range(len(d)):
                    if i!= np.argmax(d):
                        tree.remove( s[b[i]] )
                        r.update([b[i]])
                    else:
                        z=[]

        cand_list = [i for j, i in enumerate(Q) if j not in list(r)]
        #print "cand_list=", cand_list
        return cand_list



    def maxLk_interval(self, z, zs): 
        """
        print "****----in maxLk -----"
        print "z=", z,  ".....argmax=", np.argmax(z)
        print "zs=", zs, ".....argmin=", np.argmin(zs)
        print "***-....."
        """
        izmax = np.argmax(z)
        zmax  = np.max(z)

        print "izmax=", izmax
        print "zmax=", zmax


        """
        ###  (31-dec-2015):  don't remember logic of this; but it 
        ###                  is causing problems!..
        rb = [k for k in range(len(z)) if np.abs(zmax - z[k]) < 0.025  ]
        print "rb=",rb

        if len(rb)>1: 
            sb = [ zs[k] for k in rb ]
            print "sb=",sb
            izmax = np.argmin(sb)
        """

        print "final izmax=", izmax
        return izmax


#-------------------
if __name__ == '__main__':

    Q= [((6564, 6873), 0, 1.0, np.array([ 1.,  1.,  1.]), 3.0), 
        ((11668, 11974), 0, 1.0, np.array([ 1.,  1.,  1.]), 3.0), 
        ((53452, 53768), 2, 0.93600, np.array([ 0.812,  0.746,  0.936]), 2.494), 
        ((53518, 53835), 2, 0.758, np.array([ 0.758,  0.652,  0.692]), 2.10199), 
        ((53694, 54000), 2, 0.874, np.array([ 0.794,  0.874,  0.708]), 2.37599), 
        ((53724, 54039), 2, 0.89, np.array([ 0.794,  0.89 ,  0.772]), 2.456), 
        ((53733, 54039), 2, 0.87, np.array([ 0.796,  0.87 ,  0.794]), 2.46), 
        ((53810, 54115), 2, 0.93, np.array([ 0.93 ,  0.806,  0.9  ]), 2.636), 
        ((53840, 54156), 2, 0.912, np.array([ 0.912,  0.878,  0.862]), 2.652), 
        ((53849, 54156), 2, 0.906, np.array([ 0.906,  0.888,  0.88 ]), 2.6739), 
        ((53879, 54185), 2, 0.851998, np.array([ 0.838,  0.852,  0.746]), 2.43599), 
        ((53879, 54194), 2, 0.87, np.array([ 0.81 ,  0.87 ,  0.716]), 2.39599), 
        ((53888, 54194), 2, 0.87, np.array([ 0.79 ,  0.87 ,  0.748]), 2.40799), 
        ((53888, 54204), 2, 0.874, np.array([ 0.828,  0.874,  0.786]), 2.488), 
        ((53888, 54211), 2, 0.872, np.array([ 0.85 ,  0.872,  0.78 ]), 2.50198), 
        ((53918, 54233), 2, 0.743999, np.array([ 0.734,  0.71 ,  0.744]), 2.18797), 
        ((53927, 54233), 2, 0.754, np.array([ 0.702,  0.706,  0.754]), 2.16199), 
        ((53938, 54252), 2, 0.914, np.array([ 0.914,  0.816,  0.802]), 2.532), 
        ((53957, 54272), 2, 0.774, np.array([ 0.718,  0.656,  0.774]), 2.148), 
        ((53965, 54272), 2, 0.914, np.array([ 0.914,  0.912,  0.806]), 2.632), 
        ((53995, 54311), 2, 0.949996, np.array([ 0.95 ,  0.934,  0.9  ]), 2.78398), 
        ((54004, 54311), 2, 0.947995, np.array([ 0.948,  0.898,  0.874]), 2.71998), 
        ((54034, 54348), 2, 0.947995, np.array([ 0.948,  0.84 ,  0.888]), 2.676), 
        ((54043, 54348), 2, 0.937994, np.array([ 0.938,  0.826,  0.9  ]), 2.66397), 
        ((54326, 54643), 2, 0.867999, np.array([ 0.856,  0.868,  0.846]), 2.56998), 
        ((54335, 54643), 2, 0.859999, np.array([ 0.852,  0.86 ,  0.844]), 2.556), 
        ((54365, 54669), 2, 0.857998, np.array([ 0.848,  0.858,  0.828]), 2.53398), 
        ((54374, 54691), 2, 0.87, np.array([ 0.852,  0.87 ,  0.818]), 2.54), 
        ((54404, 54720), 2, 0.851998, np.array([ 0.8  ,  0.804,  0.852]), 2.456), 
        ((54404, 54727), 2, 0.853998, np.array([ 0.78 ,  0.786,  0.854]), 2.41999), 
        ((54413, 54720), 2, 0.859999, np.array([ 0.78 ,  0.754,  0.86 ]), 2.394), 
        ((54413, 54727), 2, 0.831996, np.array([ 0.752,  0.758,  0.832]), 2.342), 
        ((54959, 55263), 2, 0.723998, np.array([ 0.656,  0.586,  0.724]), 1.966), 
        ((55965, 56281), 2, 0.912, np.array([ 0.912,  0.83 ,  0.836]), 2.57798), 
        ((55974, 56281), 2, 0.926, np.array([ 0.926,  0.826,  0.834]), 2.58599), 
        ((56004, 56310), 2, 0.853998, np.array([ 0.854,  0.818,  0.76 ]), 2.43199), 
        ((56004, 56319), 2, 0.851998, np.array([ 0.852,  0.826,  0.746]), 2.42399), 
        ((56013, 56319), 2, 0.839997, np.array([ 0.84 ,  0.84 ,  0.792]), 2.472), 
        ((56052, 56358), 2, 0.845997, np.array([ 0.77 ,  0.846,  0.652]), 2.26798), 
        ((56198, 56505), 2, 0.920, np.array([ 0.92 ,  0.814,  0.824]), 2.55798), 
        ((56198, 56515), 2, 0.9240, np.array([ 0.924,  0.858,  0.83 ]), 2.612),
        ((56237, 56542), 2, 0.957996, np.array([ 0.958,  0.894,  0.932]), 2.78398), 
        ((59837, 60146), 0, 1.0, np.array([ 1.   ,  1.   ,  0.994]), 2.99398), 
        ((66519, 66828), 0, 1.0, np.array([ 1.,  1.,  1.]), 3.0), 
        ((87163, 87466), 0, 0.998, np.array([ 0.986,  0.998,  0.974]), 2.9580), 
        ((93938, 94243), 0, 0.936, np.array([ 0.9  ,  0.828,  0.936]), 2.6640), 
        ((93938, 94252), 0, 0.996, np.array([ 0.992,  0.97 ,  0.996]), 2.9580), 
        ((108611, 108920), 0, 1.0, np.array([ 0.976,  1.   ,  1.   ]), 2.976)]
 

    T = Test()

    cand_list = T.descriminate_by_exon(Q)
    

    for i in cand_list: 
        print i
