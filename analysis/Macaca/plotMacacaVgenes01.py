#!/usr/bin/env python
"""
   dnolivieri: updated 28-march-2015

     plots the loci for results.

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
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import itertools


class plotExonMap:
    def __init__(self, exonLabels, exons, xstart, xend ):
        self.height = 0.25

        self.exonLabels = exonLabels
        self.exons = exons

        self.xstart = xstart
        self.xend = xend
        self.xLower = -400
        self.xUpper =  550 #6236
        
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

            epsilon=-0.35

            ax.annotate(self.exonLabels[cnt], (x+(width)/2.,y-epsilon),
                        fontsize=10, ha='center', va='center')
            cnt+=1


        colors = 100*np.random.rand(len(patches))
        q = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.6)
        q.set_array(np.array(colors))
        ax.add_collection(q)

        
        ax.set_xlim([self.xLower, self.xUpper])
        ax.set_ylim([4, 6])



    def plot_Vexons(self, ybar, Vbar, fig, ax):
       
        x=1
        width=5
        patches=[]
        y=ybar

        cnt=0
        for p in Vbar:
            x= p[0]-xstart
            width= p[1] - p[0]
            print x, width
            rect= Rectangle( (x,y), width, self.height )
            patches.append(rect)

            epsilon=-0.35
            """
            ax.annotate(self.exonLabels[cnt], (x+(width)/2.,y-epsilon),
                        fontsize=10, ha='center', va='center')
            """
            cnt+=1


        colors = 100*np.random.rand(len(patches))
        q = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.6)
        q.set_array(np.array(colors))
        ax.add_collection(q)
        ax.set_xlim([self.xLower, self.xUpper])
        ax.set_ylim([4, 6])




    def plot_exons(self, strand, cand_exons, best_exons, fig, ax):
        patches=[]
        patchesBest = []
        y=4.0
        ymax = 4.
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

        
        ax.set_xlim([0, 6236])
        ax.set_ylim([4, 6])



# ------------------------------------------
if __name__ == '__main__':

    

    exonLabels = ['EX-4', 'EX-3', 'EX-2']

    #chr7
    S = {'TRAV8-6':(84695438,84695874),
         'TRAV40': (85118821,85119249),
         'TRAV30': (84991102,84991653),
         'TRAV25':  (84934218,84934823),
         'TRAV6':  (84496390,84496913),
         'TRAV5':  (84486442,84486950),
         'TRAV4':  (84467455,84468214),
         'TRAV27': (84971778,84972336),
         'TRAV17': (84716663,84717140),
         'TRAV19': (84732052,84732598),
         'TRAV8-4':(84564285,84564755),
         'TRAV41': (85125178,85125679)
         }



    scale = 140.
    exonLabels = [k for k in S.iterkeys() ]
    exons = [ (S[k][0]/scale, S[k][1]/scale) for k in exonLabels ]
    cand_exons = exons

    print exons
    xL = np.array([x[0] for x in exons])
    xR = np.array([x[1] for x in exons])

    print np.min(xL), np.max(xR)


    xstart= np.min(xL)
    xend=np.max(xR)
    D = plotExonMap ( exonLabels, exons, xstart, xend )


    fig, ax = plt.subplots()
    D.plot_refexons(fig, ax)


    R = {'V119RF': (376249, 376532),
         'V120RF': (401348, 401649),
         'V121RF': (488320, 488615),
         'V122RF': (500476, 500762),
         'V123RF': (519212, 519501),
         'V124RF': (529172, 529458),
         'V125RF': (565497, 565789),
         'V126RF': (578732, 579024),
         'V127RF': (591554, 591843),
         'V128RF': (596993, 597288),
         'V129RF': (617954, 618243),
         'V130RF': (640532, 640818),
         'V131RF': (652596, 652891),
         'V132RF': (666452, 666744),
         'V133RF': (672397, 672698),
         'V134RF': (690719, 691011),
         'V135RF': (728127, 728425),
         'V136RF': (739934, 740217),
         'V137RF': (749403, 749689),
         'V138RF': (761626, 761915),
         'V139RF': (764845, 765146),
         'V140RF': (800534, 800820),
         'V141RF': (810884, 811176),
         'V142RF': (839335, 839630),
         'V143RF': (872392, 872681),
         'V144RF': (880635, 880933),
         'V145RF': (891540, 891829),
         'V147RF': (955188, 955477),
         'V148RF': (959650, 959939),
         'V149RF': (967088, 967386),
         'V150RF': (978092, 978390),
         'V151RF': (987208, 987503),
         'V152RF': (1004598, 1004881),
         'V153RF': (1018152, 1018459),
         'V154RF': (1023919, 1024205),
         'V155RF': (1047022, 1047311),
         'V156RF': (1051889, 1052175),
         'V157RF': (1064813, 1065099),
         'V158RF': (1069500, 1069789),
         'V159RF': (1107467, 1107768),
         'V160RF': (1116975, 1117273),
         'V161RF': (1140579, 1140865),
         'V162RF': (1151526, 1151815),
         'V163RF': (1157941, 1158242),
         }
    

    Delta=83967455    
    exonLabels = [k for k in R.iterkeys() ]
    Vexons = [ ((R[k][0]+Delta)/scale, (R[k][1]+Delta)/scale) for k in exonLabels ]

    print Vexons
    xL = np.array([x[0] for x in exons])
    xR = np.array([x[1] for x in exons])

    print np.min(xL), np.max(xR)

    #strand = 1

    #best_exons = [0,1,2]
    #D.plot_exons(strand, cand_exons, best_exons, fig, ax)

    y=4.65
    D.plot_Vexons( y, Vexons, fig, ax)



    T={'Vs108':376251,
       'Vs109':401350,
       'Vs111':488322,
       'Vs112':500478,
       'Vs113':519213,
       'Vs114':529174,
       'Vs115':565499,
       'Vs116':578733,
       'Vs117':591555,
       'Vs118':596995,
       'Vs119':617955,
       'Vs120':640534,
       'Vs121':652598,
       'Vs122':666453,
       'Vs123':672400,
       'Vs124':690721,
       'Vs125':728129,
       'Vs126':739936,
       'Vs127':749404,
       'Vs128':761628,
       'Vs129':764847,
       'Vs130':800535,
       'Vs131':810887,
       'Vs132':839337,
       'Vs133':872393,
       'Vs134':880638,
       'Vs135':891541,
       'Vs137':955189,
       'Vs138':959651,
       'Vs139':967091,
       'Vs140':978093,
       'Vs141':987210,
       'Vs142':1004599,
       'Vs143':1023920,
       'Vs144':1047025,
       'Vs145':1051890,
       'Vs146':1064816,
       'Vs147':1069501,
       'Vs149':1107470,
       'Vs150':1116977,
       'Vs151':1140581,
       'Vs152':1157944
       }


    Delta=83967455    
    exonLabels = [k for k in T.iterkeys() ]
    Wexons = [ ((T[k]+Delta)/scale, (T[k]+300+Delta)/scale) for k in exonLabels ]

    print Wexons
    #xL = np.array([x[0] for x in exons])
    #xR = np.array([x[1] for x in exons])

    #print np.min(xL), np.max(xR)
    y=4.2
    D.plot_Vexons( y, Wexons, fig, ax)


    plt.show()

