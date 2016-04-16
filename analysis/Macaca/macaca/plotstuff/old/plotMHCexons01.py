#!/usr/bin/env python
"""
   dnolivieri  24-march-2015:
    -- taken from iface plots... 
    *** try and do the entire structure.

    future version... 
       - include annotation
       - construct directly from the fasta file.

"""


import collections
import itertools
import cPickle as pickle
import sys
import re

import numpy as np
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import itertools


class plotExonMap:
    def __init__(self, exonLabels, exons, xstart, xend ):
        self.height = 0.25

        self.exonLabels = exonLabels
        self.exons = exons

        self.xstart = xstart
        self.xend = xend

    def plot_refexons(self, fig, ax):
       
        x=1
        width=5
        patches=[]
        y=5.0

        cnt=0
        for p in self.exons:
            x= p[0]-xstart
            width= p[1] - p[0]
            print x, width
            rect= Rectangle( (x,y), width, self.height )
            patches.append(rect)

            epsilon=-0.25

            ax.annotate(self.exonLabels[cnt], (x+(width)/2., y-(self.height/2.+epsilon)), fontsize=10, ha='center', va='center')
            cnt+=1


        colors = 100*np.random.rand(len(patches))
        q = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.6)
        q.set_array(np.array(colors))
        ax.add_collection(q)


    """
    def connect_segments( s, sval, smin):


            dotted_line = plt.Line2D((x1, x2), (y1, y2), lw=1., 
                                     ls='-', marker='.', 
                                     markersize=10, 
                                     markerfacecolor='r', 
                                     markeredgecolor='r', 
                                     alpha=0.7)

    """



    def plot_exons(self, strand, cand_exons, best_exons, fig, ax):
        patches=[]
        patchesBest = []
        y=4.5
        ymax = 4.5
        height = 0.2
        xlast=0 
        epsilon=-0.35

        cnt=0
        for p in cand_exons:
            if strand<0:
                x= (self.xend - self.xstart)  - p[0] 
            else: 
                x= p[0] 

            width= strand * (p[1] - p[0])
            print x, width

            if cnt in best_exons:
                rect= Rectangle( (x,ymax), width, height )
                patchesBest.append(rect)
            else:
                rect= Rectangle( (x,y), width, height )
                patches.append(rect)


            """
            if strand<0:
                exNum = len(cand_exons)-cnt
            else:
                exNum = cnt
            """
            exNum = cnt

            if cnt in best_exons:
                ax.annotate(str(exNum), (x+(width)/2., ymax-(self.height/2.+epsilon)), fontsize=8, ha='center', va='center')
            else: 
                ax.annotate(str(exNum), (x+(width)/2., y-(self.height/2.+epsilon)), fontsize=8, ha='center', va='center')

            if np.abs(x-xlast)>100:
                y = y+0.075
            else: 
                y = y-0.075
            xlast = x
            cnt+=1

        colors = 100*np.random.rand(len(patches))
        q = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.2)
        q.set_array(np.array(colors))
        ax.add_collection(q)

        qB = PatchCollection(patchesBest, cmap=matplotlib.cm.jet, alpha=0.9)
        qB.set_array(np.array(colors))
        ax.add_collection(qB)


        ax.set_xlim([3500, 7000])
        ax.set_ylim([4, 6])



# ------------------------------------------
if __name__ == '__main__':

    
    """
    # H2-M10.6
    exons=[(36812486, 36812755), (36812986, 36813261), (36813809, 36814084)]
    xstart = 36810475
    cand_exons=[(2009, 2281), (2509, 2787), (2644, 2919),(2653, 2919),(3332, 3610), (6291, 6566), (6291, 6575),(6306, 6575), (6309, 6575)]
    """

    """
    # H2-Q7
    exons = [(35439452, 35439721),(35439909, 35440184),(35442267, 35442542)]
    cand_exons= [(1261, 1545),(2299, 2577),(2299, 2583),(2371, 2640),(2604, 2876),(3061, 3339),(7890, 8159),(7893, 8159),(8894, 9169),(8942, 9211)]

    cand_exons=[(1261, 1545),(1261, 1557),(2299, 2577),(2299, 2583),(2341, 2640),(2344, 2640),(2371, 2640),(2371, 2667),(2464, 2754),(2604, 2876),(3061, 3339),(5419, 5706),(5419, 5715),(7890, 8159),(7893, 8159),(8894, 9169),(8942, 9211)]

    xstart=35436846
    """


    exonLabels = ['EX-4', 'EX-3', 'EX-2']



    exons = [(36040478,  36040753), (36041449,  36041685), (36041868,  36042128)]
    cand_exons = [(778, 1062), (841, 1140), (844, 1140), (847, 1140), (1320, 1610), (1727, 2005), (1883, 2149), (2204, 2473), (2252, 2545), (2294, 2581), (2917, 3195), (5000, 5272), (5180, 5446)]




    xstart=36035914
    xend = 36043007


    D = plotExonMap ( exonLabels, exons, xstart, xend )


    fig, ax = plt.subplots()
    D.plot_refexons(fig, ax)

    strand = -1

    best_exons = [8]
    D.plot_exons(strand, cand_exons, best_exons, fig, ax)

    plt.show()

