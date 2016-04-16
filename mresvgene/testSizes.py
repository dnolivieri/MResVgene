#!/usr/bin/env python
"""
   dnolivieri:  23 dec 2015
     Bootstrap workflow for the MR ensemble RF code.

"""
import cPickle as pickle
import os, fnmatch
import sys


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


        
def get_existing_data(infile):
    qp = open(infile, 'rb')
    D = pickle.load(qp)
    return D

# -----------------------------------------------
if __name__ == '__main__':

    fbar= ["ighv_0.pkl", "ighv_1.pkl"]

    for ik in fbar: 
        infile='./bstrap/'+ik
        D = get_existing_data(infile)
        print len(D[2])
