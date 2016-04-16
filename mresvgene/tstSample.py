#!/usr/bin/env python
"""
   dnolivieri:  23 dec 2015
     Bootstrap workflow for the MR ensemble RF code.

"""

import collections
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import cPickle as pickle






# -----------------------------------------------
if __name__ == '__main__':

    infile="./bstrap/bckgnd_n2000.pkl"
    qp = open(infile, 'rb')
    zB = pickle.load(qp)

    jz=rnd.sample(range(0,20),3)

    print jz
    print "zB="

    zx = [zB[0][i] for i in jz]


    k=jz[0]
    print zx[0][0:4]
    print zB[0][k][0:4]

    print 
    k=jz[1]
    print zx[1][0:4]
    print zB[0][k][0:4]



    #zB[k][0:G[loci]]


    print "len(zx)=", len(zx)

    zy = zB[0][0:3]    
    print "len(zy)=", len(zy)
