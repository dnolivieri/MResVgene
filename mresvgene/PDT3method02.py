#!/usr/bin/env python
"""
   dnolivieri:  updated ...4 dec 2015
     * used to generate the "positive" multiresolution signals.
     - based upon loci.

"""
import numpy as np
import sys
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json
import cPickle as pickle
import errno
import collections

rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


class PDT3:
    def __init__(self):
        self.normI=self.normalized_AAindx()

    def normalized_AAindx(self):
        fp = open('./aaRefs/aaindex.txt','r')
        D=[]
        for lk in fp:
            q=[ float(i) for i in lk.split() ]
            D.append(q)

        Dvec=[]
        normI = []
        for j in  D:
            q= np.sum(np.array( j ))/20.
            denom=0.0
            for kp in j:
                denom+= (kp - q)*(kp - q)

            denom = np.sqrt(denom/20.)
            abar=[]
            for kp in j:
                abar.append( (kp - q)/denom )

            normI.append(abar)

        save_object(normI, r'./aaRefs/normalizedAA_Matrix.pkl')
        return normI

    def get_seq2vector(self, seq):
        Dvec={}
        cnt=0
        for q in self.normI:
            sumDseq=0.0
            for i in range(len(seq)-3):
                sumDseq+= (q[rno[seq[i]]] - q[rno[seq[i+3]]])*(q[rno[seq[i]]] - q[rno[seq[i+3]]])

            sumDseq = sumDseq/np.float(len(seq)-3)
            Dvec.update( {str(cnt): sumDseq} )
            cnt+=1
        return Dvec

    def get_lambda_seq2vector(self, seq, qLambda):
        Dvec={}
        cnt=0
        for q in self.normI:
            sumDseq=0.0
            for i in range(len(seq)- qLambda):
                pi1 = rno[seq[i]]
                pi2 = rno[seq[i+qLambda]]
                sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
            sumDseq = sumDseq/np.float(len(seq)-1)
            Dvec.update( {str(cnt): sumDseq} )
            cnt+=1
        return Dvec


    def count_sequential_doubles(self, s):
        dbar = [ s[i-1]+s[i] for i in range(len(s)) if i>0 ]
        res=collections.Counter(dbar)
        return res

    def get_freq_seq2vector(self, s):
        def merge_two_dicts(x, y):
            '''Given two dicts, merge them into a new dict as a shallow copy.'''
            z = x.copy()
            z.update(y)
            return z
        aa ='ARNDCQEGHILKMFPSTWYV'
        ## single
        dbar = [ x for x in aa] 
        S={}
        for i in dbar:
            S.update({i:0})
        T=collections.Counter(s)
        y = merge_two_dicts(S, T)        
        ybar = np.array([ y[j] for j in dbar])

        ## double.
        dbar = [ x[0]+x[1] for x in list(itertools.product(aa, repeat=2))] 
        D={}
        for i in dbar:
            D.update({i:0})
        B= self.count_sequential_doubles(s)
        z = merge_two_dicts(D, B)
        zbar = np.array([ z[j] for j in dbar ])
        
        res = np.hstack( (ybar,zbar))
        return res


## ---------------MAIN ----------------------------------
if __name__ == '__main__':


    
    P = PDT3()
