#!/usr/bin/env python
"""
   dnolivieri:  updated ...8 december 2015

      This does the prediction of Vs based 
      upon a multi-resolution training.

      - also should implement (eventually) the 
       bootstrap learning.
"""

import collections
import numpy as np
import itertools


def rhofunction(D, loci_classes):
    nlevel=2
    adapt_threshold = 0.95

    for j in range( len(D[loci_classes[0]])):  
        rho = [ D[x][j] for x in loci_classes ]
        #print j,  rho
        tau = np.array([ D[x][j].sum() for x in loci_classes ])
        #print j, tau
        loci_kmax = np.argmax(tau)
        loci_ymax = np.max(tau)
        print "j=", j, "  rho=", rho 
        print "..... tau=", tau, ".....", loci_kmax, loci_ymax
        
        mrlevel_max =np.power(2,nlevel)-1
        if ( loci_ymax  > adapt_threshold * mrlevel_max ):
                print "****found ", loci_ymax, "(>",adapt_threshold * mrlevel_max, ")....", rho[loci_kmax], rho[loci_kmax][0]

        print 




#----------------------------------
if __name__ == '__main__':



    loci_classes = ["ighv", "igkv", "iglv", "trav", "trbv", "trgv", "trdv"]

    """
    D= {'ighv': ['a1', 'a2'],
        'igkv': ['b1', 'b2'],
        'iglv': ['c1', 'c2'],
        'trav': ['d1', 'd2'],
        'trbv': ['e1', 'e2'],
        'trgv': ['f1', 'f2'],
        'trdv': ['g1', 'g2']}
    """

    D= {'trav': [np.array([ 0.91 ,  0.48 ,  0.574]), np.array([ 0.888,  0.388,  0.93 ]), 
                 np.array([ 0.918,  0.442,  0.768]), np.array([ 0.94 ,  0.434,  0.712]), 
                 np.array([ 0.924,  0.38 ,  0.522])], 
        'trgv': [np.array([ 0.496,  0.178,  0.396]), np.array([ 0.488,  0.146,  0.456]), 
                 np.array([ 0.438,  0.182,  0.294]), np.array([ 0.504,  0.214,  0.224]), 
                 np.array([ 0.422,  0.292,  0.154])], 
        'ighv': [np.array([ 0.946,  0.744,  0.672]), np.array([ 0.986,  0.986,  0.992]), 
                 np.array([ 0.986,  0.728,  0.876]), np.array([ 0.99 ,  0.706,  0.472]), 
                 np.array([ 0.956,  0.56 ,  0.27 ])], 
        'trdv': [np.array([ 0.518,  0.306,  0.538]), np.array([ 0.48 ,  0.12 ,  0.762]), 
                 np.array([ 0.438,  0.13 ,  0.49 ]), np.array([ 0.356,  0.126,  0.334]), 
                 np.array([ 0.292,  0.146,  0.194])], 
        'iglv': [np.array([ 0.74 ,  0.198,  0.542]), np.array([ 0.736,  0.118,  0.752]), 
                 np.array([ 0.746,  0.208,  0.568]), np.array([ 0.802,  0.158,  0.52 ]), 
                 np.array([ 0.844,  0.35 ,  0.498])], 
        'trbv': [np.array([ 0.696,  0.274,  0.544]), np.array([ 0.682,  0.236,  0.648]), 
                 np.array([ 0.704,  0.208,  0.748]), np.array([ 0.762,  0.258,  0.786]), 
                 np.array([ 0.828,  0.29 ,  0.62 ])], 
        'igkv': [np.array([ 0.792,  0.216,  0.37 ]), np.array([ 0.774,  0.204,  0.534]), 
                 np.array([ 0.75 ,  0.38 ,  0.324]), np.array([ 0.818,  0.368,  0.298]), 
                 np.array([ 0.74 ,  0.378,  0.266])]}

    rhofunction(D, loci_classes)
