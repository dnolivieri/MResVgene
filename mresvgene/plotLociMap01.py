#!/usr/bin/env python
"""
   dnolivieri  11-enero-2016
     - For making maps of the loci.
"""
import collections
import itertools
import sys
import re
import cPickle as pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
import matplotlib as mpl
import seaborn as sns
from numpy.random import randn
import scipy
import numpy as np
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import itertools


class plotExonMap:
    def __init__(self, xstart, xend ):
        self.height = 0.25
        self.xstart = xstart
        self.xend = xend


    def get_references(self, inFile):
        for record in SeqIO.parse(inFile, "fasta"):
            rec_name=str(record.name)
            rec_desc=str(record.description)
            print rec_desc


    def plot_exons(self, cand_exons, fig, ax):
        patches=[]
        patchesA = []


        y=4.5
        ymax = 4.5
        height = 0.2
        xlast=0 
        epsilon=-0.35

        cnt=0
        for p in cand_exons:
            x= p[0] 
            width= (p[1] - p[0])
            rect= Rectangle( (x,y), width, height )
            patches.append(rect)
            delta=-0.2

            #ax.annotate("ex1", (x+(width)/2., y-(self.height/2.+delta)), fontsize=10, ha='center', va='center')
            rect= Rectangle( (x,y), width, height, color='blue',  alpha=0.9)
            patchesA.append(rect)

            #ax.annotate(str(exNum), (x+(width)/2., y-(self.height/2.+epsilon)), fontsize=8, ha='center', va='center')
            y = y-0.075
            
            xlast = x
            cnt+=1


        colorsA = 100 * np.ones(len(patchesA), dtype=np.int)
        qA = PatchCollection(patchesA, cmap=matplotlib.cm.jet, alpha=0.6)
        qA.set_clim([5,50])
        qA.set_array(np.array(colorsA))
        ax.add_collection(qA)


        ax.set_xlim([xstart, xend])
        ax.set_ylim([0, 6])



# ------------------------------------------
if __name__ == '__main__':

    infile="./IMGTrefs/Macaca_mulata_IMGT.fasta"

    xstart=53500
    xend = 67000

    D = plotExonMap ( xstart, xend )
    D.get_references(infile)


    """
    cand_exons = [(6564, 6873), (11668, 11974), (53452, 53768)]
    cand_exons=[(6564, 6873), (11668, 11974), (53995, 54311)]

    xstart=53500
    xend = 67000

    D = plotExonMap ( xstart, xend )

    fig, ax = plt.subplots()
    D.plot_exons( cand_exons, fig, ax)

    plt.show()
    """
