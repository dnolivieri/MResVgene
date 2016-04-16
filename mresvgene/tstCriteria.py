#!/usr/bin/env python
"""
  dnolivieri:  updated ...7-jan-2016
   -- results from the analysis of all 
      sequences that should be included.
"""
import numpy as np
import re
import time
import os, fnmatch
import sys
import itertools

nlevel=2


def MRscore(pm):
    p = {}
    cnt=0
    for k in range(nlevel):
        p.update({k:[]})
        for m in range(np.power(2,k)):
            p[k].append(pm[cnt])
            cnt+=1
    sum_score = pm.sum()
    sigma=0.25
    sbar=0.
    for k in range(1,nlevel):
        for m in range(np.power(2,k)):
            Delta = np.abs(p[0][0] - p[k][m])
            sdelta = 1.0- np.exp( -np.power(Delta,2)/sigma ) 
            sbar+= sdelta
            cnt+=1

    tot_score = sum_score - sbar 
    return tot_score

def get_result(q):
    threshold = 0.725
    N_mrlevel = np.power(2,nlevel)-1
    qMR = np.array([ MRscore(q[i])  for i  in range(N_mrlevel) ])
    iqMRmax = np.argmax(qMR)
    qMRmax  = np.max(qMR)
    
    qth = threshold * N_mrlevel

    #if (q[iqMRmax].sum() > qth)
    print qMR,qMRmax, qth,  q[iqMRmax].sum(), "|", q[3],  
    if (qMRmax > qth):
        print "***",
        sbar=True
        rx=q[iqMRmax]
        r_ix=np.where(np.abs(rx-rx[0])>0.15)[0]
        if len(r_ix) >= int( (N_mrlevel-1)/2.):
            for k in range(len(r_ix)): 
                print "|", r_ix,  rx[ r_ix[k] ], "|", rx, 
                if  rx[ r_ix[k] ] < 0.5:
                    sbar=False

        print sbar
    else:

        sbar=False
        print sbar
        """
        if sbar==True:
            seqs.append( ( seqpos[j],seqbuffer[j], loci_kmax, np.max(rho[loci_kmax]), Pr_outclass, np.sum(Pr_outclass) ) )
        """

def run(Q): 

    for q in Q:
        get_result(q)
        



#-------------------
if __name__ == '__main__':
    """
    1==yes
    0==no
    -1 == not sure
    """
    Q = [[np.array([ 0.201,  0.014, 0.495]),  np.array([ 0.46 , 0.847, 0.675]), np.array([ 0.949,  0.964, 0.959]), 1],
         [np.array([ 0.201,  0.014, 0.495]),  np.array([ 0.46 , 0.847, 0.675]), np.array([ 0.949,  0.664, 0.959]), 1],
         [np.array([ 0.213,  0.011, 0.205]),  np.array([ 0.4  , 0.263, 0.261]), np.array([ 0.988,  0.991, 0.834]), 1],
         [np.array([ 0.094,  0.011, 0.012]),  np.array([ 0.222, 0.297, 0.112]), np.array([ 0.888,  0.93 , 0.175]), 0],
         [np.array([ 0.158,  0.006, 0.165]),  np.array([ 0.274, 0.209, 0.197]), np.array([ 0.962,  0.929, 0.636]), 1], 
         [np.array([ 0.219,  0.009, 0.018]),  np.array([ 0.365, 0.256, 0.083]), np.array([ 0.816,  0.882, 0.234]), 0],
         [np.array([ 0.291,  0.004,  0.14]),  np.array([ 0.397, 0.178, 0.197]), np.array([ 0.803,  0.825, 0.652]), 1],
         [np.array([ 0.118,  0.022,  0.01]),  np.array([ 0.929, 0.751, 0.523]), np.array([ 0.664,  0.217, 0.036]),-1],
         [np.array([ 0.597,  0.859, 0.305]),  np.array([ 0.323, 0.154, 0.028]), np.array([ 0.043,  0.147, 0.055]), 0],
         [np.array([ 0.677,  0.948, 0.428]),  np.array([ 0.291, 0.129, 0.022]), np.array([ 0.074,  0.1  , 0.034]),-1],
         [np.array([ 1.,     0.996, 0.848]),  np.array([ 0.559, 0.15 , 0.213]), np.array([ 0.576,  0.17 , 0.218]), 1],
         [np.array([ 1.,     0.658,   1. ]),  np.array([ 0.808, 0.188, 0.057]), np.array([ 0.555,  0.141, 0.161]), 1],
         [np.array([ 0.998,  0.851, 0.997]),  np.array([ 0.715, 0.038, 0.087]), np.array([ 0.51 ,  0.033, 0.128]), 1],
         [np.array([ 0.978,  0.751,    1.]),  np.array([ 0.388, 0.038, 0.047]), np.array([ 0.262,  0.018, 0.173]), 1],
         [np.array([ 0.801,  0.026, 0.502]),  np.array([ 0.997, 0.788, 0.984]), np.array([ 0.867,  0.145, 0.781]), 1] ]
    run(Q)



