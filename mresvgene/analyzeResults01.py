#!/usr/bin/env python
"""
   dnolivieri:  04 jan 2016
      analysis of the results from the MR code....

"""

import collections
import numpy as np
import matplotlib.pyplot as plt
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

from scipy import *
import struct
import re
import json
import cPickle as pickle
from collections import defaultdict
from copy import deepcopy

import timeit
import operator

from time import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection)




from PDT3method import PDT3
from MResStruct02 import MResStruct

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



class AnalyzeMRres:
    def __init__(self, S, loci_classes):
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




    def get_sequences(self):
        """
        X = ...
        y = ....
        n_samples, n_features = X.shape
        n_neighbors = 30
        """
        pass

    
    def plotx(self):
        print("Computing t-SNE embedding")
        tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
        t0 = time()
        X_tsne = tsne.fit_transform(X)
        plot_embedding(X_tsne,"t-SNE embedding of the digits (time %.2fs)" % (time() - t0))
        plt.show()

        
        def plot_embedding(X, title=None):
            # Scale and visualize the embedding vectors
            x_min, x_max = np.min(X, 0), np.max(X, 0)
            X = (X - x_min) / (x_max - x_min)
            
            plt.figure()
            ax = plt.subplot(111)
            for i in range(X.shape[0]):
                plt.text(X[i, 0], X[i, 1], str(digits.target[i]),
                         color=plt.cm.Set1(y[i] / 10.),
                         fontdict={'weight': 'bold', 'size': 9})

            if hasattr(offsetbox, 'AnnotationBbox'):
                # only print thumbnails with matplotlib > 1.0
                shown_images = np.array([[1., 1.]])  # just something big
                for i in range(digits.data.shape[0]):
                    dist = np.sum((X[i] - shown_images) ** 2, 1)
                    if np.min(dist) < 4e-3:
                        # don't show points that are too close
                        continue
                    shown_images = np.r_[shown_images, [X[i]]]
                    imagebox = offsetbox.AnnotationBbox(
                        offsetbox.OffsetImage(digits.images[i], cmap=plt.cm.gray_r),
                        X[i])
                    ax.add_artist(imagebox)
            plt.xticks([]), plt.yticks([])
            if title is not None:
                plt.title(title)
                        


    def run(self):
        print "forming positive npz files"
        self.form_featurevecs()




# -----------------------------------------------
if __name__ == '__main__':


    Vs_wgs_prediction = 'wgs_primates.json'
    json_data2=open( Vs_wgs_prediction )
    P = json.load(json_data2)
    json_data2.close()

