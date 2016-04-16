#!/usr/bin/env python
"""
   dnolivieri:  updated ...24 dec 2015
      - this is used in the bootstrap code.
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
import collections

from PDT3method import PDT3
from MResStruct02 import MResStruct


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class getMRFeatureVectors:
    def __init__(self,S):
        self.S = S
        self.pdt3 = PDT3()
        self.nlevel = 2
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
                P  = self.Mstruct.get_pyramid(seq_data)
                knt=0
                for k in range(len(P)):
                    qLambda = self.nlevel - k
                    #qLambda = 1  ## if the different levels are not multires.
                    for kseq in P[k]:
                        T = self.pdt3.get_lambda_seq2vector( kseq, qLambda )
                        Tx = [ T[str(x)]  for x in range(len(T)) ]
                        D[knt].append(Tx)
                        print cnt, k, knt, " ....", Tx[0:3],   kseq
                        knt+=1
                print
                cnt+=1
                if cnt>1e9:
                    break

        return D



    def background_from_fasta(self, infile, maxseqs):
        cnt=0
        D = {}
        for k in range(np.power(2,self.nlevel)-1):
            D.update( {k:[]} )
        for record in SeqIO.parse(infile, "fasta"):
            descObject=PyPro.GetProDes(record.seq.tostring())
            seq_data = record.seq.tostring()
            if ('X' not in seq_data ) and ('Z' not in seq_data) and ('B' not in seq_data):
                P  = self.Mstruct.get_pyramid(seq_data)
                knt=0
                for k in range(len(P)):
                    qLambda = self.nlevel - k
                    #qLambda = 1  ## if the different levels are not multires.
                    for kseq in P[k]:
                        T = self.pdt3.get_lambda_seq2vector( kseq, qLambda )
                        Tx = [ T[str(x)]  for x in range(len(T)) ]
                        D[knt].append(Tx)
                        print cnt, k, knt, " ....", Tx[0:3],   kseq
                        knt+=1
                print
                cnt+=1
                if (cnt>maxseqs) or (cnt>1e9):
                    break

        return D





    def hybrid_descriptors_from_fasta(self, infile, maxseqs=100000):
        cnt=0
        D = {}
        for k in range(np.power(2,self.nlevel)-1):
            D.update( {k:[]} )
        for record in SeqIO.parse(infile, "fasta"):
            descObject=PyPro.GetProDes(record.seq.tostring())
            seq_data = record.seq.tostring()
            if ('X' not in seq_data ) and ('Z' not in seq_data) and ('B' not in seq_data):
                P  = self.Mstruct.get_pyramid(seq_data)
                knt=0
                for k in range(len(P)):
                    qLambda = self.nlevel - k
                    #qLambda = 1  ## if the different levels are not multires.
                    for kseq in P[k]:
                        T={}
                        if k==0: 
                            T = self.pdt3.get_freq_seq2vector(kseq)
                        else:
                            #T = self.pdt3.get_lambda_seq2vector( kseq, qLambda )
                            T = self.pdt3.get_freq_seq2vector(kseq)

                        Tx = [ T[str(x)]  for x in range(len(T)) ]
                        D[knt].append(Tx)
                        print cnt, k, knt, " ....", Tx[0:6],   kseq
                        knt+=1
                print
                cnt+=1
                if (cnt>maxseqs) or (cnt>1e9):
                    break

        return D





## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    
    Vs_Loci = 'Vs_Loci_fasta.json'
    json_data=open( Vs_Loci )
    S = json.load(json_data)
    json_data.close()





