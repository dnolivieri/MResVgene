#!/usr/bin/env python
"""
   dnolivieri:  updated ...29-jan-2016
     - test the occurences...
"""
import collections
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


class Test:
    def __init__(self, inFile):
        self.inFile= inFile


    def get_occurrence_stats(self):
        for record in SeqIO.parse(self.inFile, "fasta"):
            seq= record.seq
            #res= self.count_singles(seq)
            res= self.count_sequential_doubles(seq)
            print record.description, "......", res.most_common(3), res.most_common(3)[0][1]

    def count_singles(self, s):
        res=collections.Counter(s)
        #print res
        #print res.most_common(2), res.most_common(2)[0][1]
        return res


    def count_sequential_doubles(self, s):
        dbar = [ s[i-1]+s[i] for i in range(len(s)) if i>0 ]
        res=collections.Counter(dbar)
        #res.most_common(2)[0][1]
        #print res
        return res



    def form_feature_vec(self, s):
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
        T=self.count_singles(s)
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


    def run(self):
        self.get_occurrence_stats()
        """
        self.count_singles(s)
        self.count_sequential_doubles(s)
        """
#-------------------
if __name__ == '__main__':

    infile = "./bstrap/trdv_1.fasta"
    T = Test(infile)
    #T.run()

    s="NLITQLVTSVTKKKGNTAFLECQIKSSTFKKNVCIHWYHQKPDQPLKRILYISSNENVVH"
    r = T.form_feature_vec(s)
    
    print r
    print len(r)
