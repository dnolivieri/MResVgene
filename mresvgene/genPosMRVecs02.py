#!/usr/bin/env python
"""
   dnolivieri:  updated ...4 dec 2015
     * used to generate the "positive" multiresolution signals.
     - based upon loci.

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
import DRmethods01 as dr

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)




class genPosMRFeatureVec:
    def __init__(self,S, loci_classes):
        self.S = S
        self.desc_method='PDT'
        self.loci_classes = loci_classes
        self.pdt3 = PDT3()
        self.DR = dr.DRmethods(100)
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
                S  = self.Mstruct.get_pyramid(seq_data)
                knt=0
                for k in range(len(S)):
                    qLambda = self.nlevel - k
                    #qLambda = 1   ## the other doesn't work...
                    for kseq in S[k]:
                        #T = self.pdt3.get_seq2vector( kseq )
                        T = self.pdt3.get_lambda_seq2vector( kseq, qLambda )
                        Tx = [ T[x]  for x in T.iterkeys() ]
                        D[knt].append(Tx)
                        print cnt, k, knt, " ....", Tx[0:3],   kseq
                        knt+=1
                print
                cnt+=1
                if cnt>1e9:
                    break

        return D

    def form_featurevecs(self):
        for loci in self.loci_classes:
            loci_file = self.S['All'][loci]
            Dbar=self.descriptors_from_fasta(loci_file)

            #Dx = self.reduce_all_MR_dims( Dbar)

            pkloutfile=loci_file.replace(".fasta", ".pkl")
            print pkloutfile
            save_object(Dbar, pkloutfile)



    def reduce_all_MR_dims(self, Dbar):
        ## PUT HERE:  dim-reduction technique.
        print type(Dbar), len(Dbar)
        Dx = {}
        for k in Dbar.iterkeys():
            Dx.update({k:[]})

        for k in Dbar.iterkeys():
            X = np.array(Dbar[k])
            Xnew = self.DR.apply_DR(X)
            print type(X), X.shape, Xnew.shape

        return Dx


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
    loci_classes=[ 'ighv', 'iglv', 'igkv', 'trav','trbv','trgv', 'trdv']

    #loci_classes=[ 'ighv', 'igkv' ]
    loci_classes=[ 'bkg1', 'bkg2' ]
    #loci_classes=[ 'bkg3']

    P = genPosMRFeatureVec( S, loci_classes)
    P.run()




