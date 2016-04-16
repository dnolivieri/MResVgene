#!/usr/bin/env python
"""
   dnolivieri  09-may-2015:

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
    def __init__(self, xstart, xend ):
        self.height = 0.25
        self.xstart = xstart
        self.xend = xend



    def plot_exons(self, strand, cand_exons, class_exons, best_exons, fig, ax):
        patches=[]
        patchesBest = []

        patchesA = []
        patchesB = []
        patchesC = []


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


            delta=-0.2
            if  class_exons[cnt]==0: 
                ax.annotate("ex1", (x+(width)/2., y-(self.height/2.+delta)), fontsize=10, ha='center', va='center')
                rect= Rectangle( (x,y), width, height, color='blue',  alpha=0.9)
                patchesA.append(rect)
            elif class_exons[cnt]==1: 
                ax.annotate("ex2", (x+(width)/2., y-(self.height/2.+delta)), fontsize=10, ha='center', va='center')
                rect= Rectangle( (x,y), width, height, color='red',  alpha=0.9)
                patchesB.append(rect)
            else: 
                ax.annotate("ex3", (x+(width)/2., y-(self.height/2.+delta)), fontsize=10, ha='center', va='center')
                rect= Rectangle( (x,y), width, height, color='green',  alpha=0.9)
                patchesC.append(rect)

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


        """
        colors = 100*np.random.rand(len(patches))
        q = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.2)
        q.set_array(np.array(colors))
        ax.add_collection(q)

        qB = PatchCollection(patchesBest, cmap=matplotlib.cm.jet, alpha=0.9)
        qB.set_array(np.array(colors))
        ax.add_collection(qB)
        """

        colorsA = 10 * np.ones(len(patchesA), dtype=np.int)
        print "colorsA=", colorsA
        qA = PatchCollection(patchesA, cmap=matplotlib.cm.jet, alpha=0.6)
        qA.set_clim([5,50])
        qA.set_array(np.array(colorsA))
        ax.add_collection(qA)


        colorsB = 30 * np.ones(len(patchesB), dtype=np.int)
        print "colorsB=", colorsB
        qB = PatchCollection(patchesB, cmap=matplotlib.cm.jet, alpha=0.6)
        qB.set_clim([5,50])
        qB.set_array(np.array(colorsB))
        ax.add_collection(qB)

        colorsC = 49 * np.ones(len(patchesC), dtype=np.int)
        print "colorsC=", colorsC
        qC = PatchCollection(patchesC, cmap=matplotlib.cm.jet, alpha=0.6)
        qC.set_clim([5,50])
        qC.set_array(np.array(colorsC))
        ax.add_collection(qC)

        ax.set_xlim([xstart, xend])
        ax.set_ylim([4, 6])



# ------------------------------------------
if __name__ == '__main__':

    

    cand_exons = [(778, 1062), (841, 1140), (844, 1140), (847, 1140), (1320, 1610), (1727, 2005), (1883, 2149), (2204, 2473), (2252, 2545), (2294, 2581), (2917, 3195), (5000, 5272), (5180, 5446)]

    cand_exons= [(15475, 15756), (15486, 15773), (15547, 15819), (16016, 16294), (16870, 17148), (16879, 17148), (27842, 28111), (28260, 28550), (28263, 28550), (28275, 28574), (28284, 28574), (28299, 28574), (46016, 46285), (46132, 46422), (46171, 46440), (46174, 46440), (46201, 46494), (46340, 46612), (46615, 46884), (46814, 47092), (47684, 47962), (47693, 47962), (60304, 60573)]


    cand_exons= [(15486, 15773), (16016, 16294), (16870, 17148), (27842, 28111), (28275, 28574), (46016, 46285), (46340, 46612), (46615, 46884), (46814, 47092), (47684, 47962), (60304, 60573)]






    #step 0
    cand_exons= [(15475, 15756), (15486, 15773), (15547, 15819), (16016, 16294), (16870, 17148), (16879, 17148), (27842, 28111), (28260, 28550), (28263, 28550), (28275, 28574), (28284, 28574), (28299, 28574), (46016, 46285), (46132, 46422), (46171, 46440), (46174, 46440), (46201, 46494), (46340, 46612), (46615, 46884), (46814, 47092), (47684, 47962), (47693, 47962), (60304, 60573)]
    class_exons= [0, 2, 0, 1, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 1, 2, 2, 2]


    # step 1

    cand_exons= [(15486, 15773), (16016, 16294), (16870, 17148), (27842, 28111), (28275, 28574), (46016, 46285), (46340, 46612), (46615, 46884), (46814, 47092), (47684, 47962), (60304, 60573)]
    class_exons= [2, 1, 2, 0, 2, 2, 0, 2, 1, 2, 2]

    """

    # step 2
    cand_exons= [(15486, 15773), (16016, 16294), (16870, 17148), (27842, 28111), (46016, 46285), (46340, 46612), (46615, 46884), (46814, 47092), (47684, 47962), (60304, 60573)]
    class_exons= [2, 1, 2, 0, 2, 0, 2, 1, 2, 2]
    """

    # step 3

    cand_exons= [(15486, 15773), (16016, 16294), (16870, 17148), (27842, 28111), (46016, 46285), (46340, 46612), (46814, 47092), (47684, 47962), (60304, 60573)]
    class_exons= [2, 1, 2, 0, 2, 0, 1, 2, 2]
    



    xstart=45400
    xend = 48500

    D = plotExonMap ( xstart, xend )

    fig, ax = plt.subplots()
    strand=1

    best_exons= [100]

    D.plot_exons(strand, cand_exons, class_exons, best_exons, fig, ax)

    plt.show()

