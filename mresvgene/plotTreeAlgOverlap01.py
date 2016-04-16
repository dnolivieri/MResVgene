#!/usr/bin/env python
"""
   dnolivieri  5-enero-2016
     - modified from plotTreeOverlap.
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

    

    cand_exons = [(6564, 6873), (11668, 11974), (53452, 53768), (53518, 53835), (53694, 54000), (53724, 54039), (53733, 54039), 
                  (53810, 54115), (53840, 54156), (53849, 54156), (53879, 54185), (53879, 54194), (53888, 54194), (53888, 54204), (53888, 54211), 
                  (53918, 54233), (53927, 54233), (53938, 54252), (53957, 54272), (53965, 54272), (53995, 54311), (54004, 54311), (54034, 54348), 
                  (54043, 54348), (54326, 54643), (54335, 54643), (54365, 54669), (54374, 54691), (54404, 54720), (54404, 54727), (54413, 54720), 
                  (54413, 54727), (54959, 55263), (55965, 56281), (55974, 56281), (56004, 56310), (56004, 56319), (56013, 56319), (56052, 56358), 
                  (56198, 56505), (56198, 56515), (56237, 56542), (59837, 60146), (66519, 66828), (87163, 87466), (93938, 94243), (93938, 94252), (108611, 108920)]   



    cand_exons=[(6564, 6873), (11668, 11974), (53995, 54311), (54374, 54691), (54959, 55263), (56237, 56542), (59837, 60146), (66519, 66828), (87163, 87466), (93938, 94252), (108611, 108920)]

    xstart=53500
    xend = 67000

    D = plotExonMap ( xstart, xend )

    fig, ax = plt.subplots()
    D.plot_exons( cand_exons, fig, ax)

    plt.show()

