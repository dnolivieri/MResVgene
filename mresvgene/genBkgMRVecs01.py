#!/usr/bin/env python
"""
   dnolivieri:  updated ...7 dec 2015
    * used to generate the multiresolution background signals
"""
import pylab as pl
import numpy as np
import sys
from numpy import genfromtxt, savetxt
import time
import timeit
import itertools
import cPickle as pickle
from Bio import SeqIO
from Bio.Seq import Seq
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json

from PDT3method import PDT3
from MResStruct02 import MResStruct

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


class genBkgMRFeatureVec:
    def __init__(self,S):
        self.S = S
        self.desc_method='PDT'
        self.pdt3 = PDT3()
        self.nlevel = 3
        self.Mstruct = MResStruct(self.nlevel)


    def descriptors_from_fasta(self, infile):
        cnt=0
        D = {}
        for k in range(np.power(2,self.nlevel)-1):
            D.update( {k:[]} )
        for record in SeqIO.parse(infile, "fasta"):
            descObject=PyPro.GetProDes(record.seq.tostring())
            seq_data = record.seq.tostring()
            if ('X' not in seq_data ) and ('Z' not in seq_data) and ('B' not in seq_data):
                S  = self.Mstruct.get_pyramid(seq_data)
                knt=0
                for k in range(len(S)):
                    for kseq in S[k]:
                        T = self.pdt3.get_seq2vector( kseq )
                        Tx = [ T[x]  for x in T.iterkeys() ]
                        D[knt].append(Tx)
                        print cnt, k, knt, " ....", Tx[0:4],   kseq
                        knt+=1
                print
                cnt+=1
                if cnt>1e9:
                    break
        return D

    def form_featurevecs(self):
        bkgnd_file = './dataVgeneDB/bkg100.fasta'
        Dbar=self.descriptors_from_fasta(bkgnd_file)
        npzoutfile= bkgnd_file.replace(".fasta", ".pkl")
        print npzoutfile
        save_object(Dbar, npzoutfile)

        """
        bkgnd_file = './dataVgeneDB/bkg.fasta'
        Dbar=self.descriptors_from_fasta(loci_file)
        npzoutfile=loci_file.replace(".fasta", ".pkl")
        print npzoutfile
        save_object(Dbar, npzoutfile)
        """

    def run(self):
        print "forming positive npz files"
        self.form_featurevecs()

        """
        ## this is the correct way to read a locus file.
        qp = open('./dataVgeneDB/ighv.pkl', 'rb')
        D = pickle.load(qp)
        print D[1]
        """

## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    
    Vs_Loci = 'Vs_Loci_fasta.json'
    json_data=open( Vs_Loci )
    S = json.load(json_data)
    json_data.close()


    B = genBkgMRFeatureVec( S )
    B.run()


