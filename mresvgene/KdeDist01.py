#!/usr/bin/env python
"""
   dnolivieri: (11jan2016):
    - this does the iteration selection step for next iteration.
 """

import shutil
import numpy as np
import time
import bisect 
import os, fnmatch
import glob
import sys
import re

import itertools
import cPickle as pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numpy.random import randn
import scipy


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



class KDEProbDistribution: 
    def __init__(self, D, loci_classes):
        self.D = D
        self.loci_classes = loci_classes


    def get_kde_struct(self, show_plot=False):
        def get_median( pbar, indx ): 
            x,y = pbar.get_lines()[indx].get_data()
            cdf = scipy.integrate.cumtrapz(y, x, initial=0)
            nearest_05 = np.abs(cdf-0.5).argmin()
            x_median = x[nearest_05]
            y_median = y[nearest_05]
            return x_median, y_median

        X_Medians=[]

        mpl.rc("figure", figsize=(8, 4))
        colors = ['r', 'y', 'g','g','g','g', 'g']
        cnt=0
        pbar=[]
        for loci in self.loci_classes:
            p = [ x[0] for x in self.D[loci] ]
            if len(p)==0: 
                pr1=np.ones(2)
            else:
                pr1 = np.array(p)

            sns.set_palette("hls")
            sns.set_style("whitegrid")
            #sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.0})
            sns.set_context("notebook", font_scale=1.2, rc={"lines.linewidth": 2.6})

            p1=sns.kdeplot(pr1, label=loci);
            x_median, y_median = get_median( p1, cnt )
            #print x_median, y_median 
            X_Medians.append(x_median)
            plt.vlines(x_median, 0, y_median, colors[cnt], linestyle='--')
            pbar.append(p1)
            cnt+=1

        if show_plot:
            plt.xlim([1.0,3.5])

            #plt.ylim([0,6.0])
            plt.ylim([0,6])
            plt.legend()
            plt.show()
            
        return X_Medians




# ------------------------------------------
if __name__ == '__main__': 
    K = KDEProbDistribution(0)






    """
    p2=sns.kdeplot(pr2)
    x2_median, y2_median = get_median( p2, 1 )
    print x2_median, y2_median
    plt.vlines(x2_median, 0, y2_median, 'y', linestyle='--')
    X_Medians.append(x2_median)
    
    p3=sns.kdeplot(pr3)
    x3_median, y3_median = get_median( p3, 2 )
    print  x3_median, y3_median
    plt.vlines(x3_median, 0, y3_median, 'g', linestyle='--')
    X_Medians.append(x3_median)
    

    ## saves the plot that can be reloaded.
    >>> q=open('kdedist.pkl', 'rb')
    >>> p=pickle.load(q)

    outfile="./bstrap/bstrp_kdedist"+str(self.cntIter) +".pkl"
    save_object(p3, outfile)
    """
