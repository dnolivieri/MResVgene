#!/usr/bin/env python
"""
   dnolivieri  24-march-2015:
    -- taken from iface plots... 
    *** try and do the entire structure.


Instead of hiding each element, you can hide the whole axis:

frame1.axes.get_xaxis().set_visible(False)
frame1.axes.get_yaxis().set_visible(False)
Or, you can set the ticks to an empty list:

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])

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

import MHCmodel02 as mhcMod
import TCRmodel01 as tcrMod
import mhcComplex01 as trMH

fig, ax = plt.subplots()

class ppPairData: 
    def __init__(self, pair, trans_pep1, trans_pep2,  p_tcr, p_mhc, trans_qcontacts  ):
        self.pair= pair
        self.trans_pep1 = trans_pep1
        self.trans_pep2 = trans_pep2
        self.trans_qcontacts = trans_qcontacts
        self.p_tcr = p_tcr
        self.p_mhc = p_mhc


class plotPeptideMap:
    def __init__(self, work_dir, bStruct ):
        self.height = 0.5
        self.bStruct = bStruct
        self.TRAV_name="TRAV"
        self.MH_name="MHC"       

        self.work_dir = work_dir

        print self.bStruct.pairs
        self.plot_map()


    def plot_map(self):
        
        yH=[9.0, 4.]
        for kindx in range(len(self.bStruct.pairs)):
            print "***", kindx
            tr_name, mh_name, tr_pep, mh_pep, qcontacts, qvals = self.extract_data(kindx)
            if kindx==0:
                self.plot_segments(yH[kindx], tr_name, mh_name, tr_pep, mh_pep, qcontacts, qvals)
            else:
                qcontacts_rev = [ [v,k] for k,v in qcontacts ]
                self.plot_segments(yH[kindx], mh_name, tr_name, mh_pep, tr_pep, qcontacts_rev, qvals)



        #plt.subplots_adjust(left=0.05, right=0.96, bottom=0.06, top=0.96)
        #plt.show()



        
        #frame1 = plt.gca()

        frame1 = fig.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        
        fig.set_size_inches(18.5,6.5)
        fig.subplots_adjust(left=0.02, right=0.97, bottom=0.06, top=0.96)
        #fig.savefig(self.work_dir + "/pepMap_" + self.work_dir.split("_")[1] +".png",dpi=100)


    def extract_data(self, kindx):
        tr=[ x for x in self.bStruct.TR_peptides[kindx] if [] not in x]
        mh=[ x for x in self.bStruct.MH_peptides[kindx] if [] not in x]
        tr_pep=[]
        for x in tr: 
            try:
                tr_pep.append(list(itertools.chain(*x)))
            except: 
                tr_pep.append(x)
        mh_pep=[]
        for x in mh: 
            try:
                mh_pep.append(list(itertools.chain(*x)))
            except: 
                mh_pep.append(x)
        qcontacts = [  [k,v] for k,v,l in self.bStruct.contacts[kindx]]
        qvals = [  l for k,v,l in self.bStruct.contacts[kindx]]


        tr_name = self.bStruct.p_tcr[kindx][0]
        mh_name = self.bStruct.p_mhc[kindx][0]

        
        print [x for x in self.bStruct.p_tcr[kindx]]
        print self.bStruct.p_tcr[kindx][0]


        cdr_filter=True

        if cdr_filter:
            s = [k[2] for k in self.bStruct.p_tcr[kindx]  ]
            s_mask = [re.match(r'CDR[1-2]', x)!=None for x in s]
            ##[(i,v) for i,v in enumerate(p) if v==True]

            sindx=[i for i,v in enumerate(s_mask) if v==True]
            qcontacts = [ qcontacts[j] for j in sindx]
            qvals = [ qvals[j] for j in sindx]
            
            #print "s=", s
            #print "s_mask=", s_mask
            #print "sindx=", sindx
            print "qcontacts=", qcontacts

            r = [k[2] for k in self.bStruct.p_mhc[kindx]  ]
            sbar = [s[j] for j in sindx]
            rbar = [r[j] for j in sindx]
            #print "sbar,rbar =", sbar, rbar 

        return tr_name, mh_name, tr_pep, mh_pep, qcontacts, qvals


    def plot_segments(self, ybar, pep1_name, pep2_name, peptides_1, peptides_2, qcontacts, qvals):
        x=1
        width=5
        patches=[]
        y=ybar
        #ax.annotate(pep1_name, (10,y-0.5), fontsize=10, ha='center', va='center')

        for p in peptides_1: 
            x= p[0]
            width= p[1] - p[0]
            rect= Rectangle( (x,y), width, self.height )
            patches.append(rect)
            ax.annotate(p[0], (x,y-0.1), fontsize=8, ha='center', va='center')
            ax.annotate(p[1], (x+width,y-0.1), fontsize=8, ha='center', va='center')

        y=y-3.5
        #ax.annotate(pep2_name, (10,y+0.9), fontsize=10, ha='center', va='center')

        x_last=0

        for q in peptides_2: 
            x= q[0]
            width=q[1] - q[0]
            rect= Rectangle( (x,y), width, self.height )
            patches.append(rect)
            if (x-x_last < 4):
                ax.annotate(q[0], (x,y+0.1+self.height), fontsize=8, ha='center', va='center')
                ax.annotate(q[1], (x+width,y+0.1+self.height), fontsize=8, ha='center', va='center')
            else: 
                ax.annotate(q[0], (x,y+0.1+self.height), fontsize=8, ha='center', va='center')
                ax.annotate(q[1], (x+width,y+0.1+self.height), fontsize=8, ha='center', va='center')
                
            x_last = x+width


            
        def connect_segments( s, sval, smin):
            x1, y1 = s[0], ybar+.25
            x2, y2 = s[1], y1-3.5

            r1= x1+(x2-x1)/2.
            r2= y2+ np.abs(y2-y1)/2.
            #print r1,r2
            #print "sval=",sval,  '{:.1f}'.format(sval)
            ax.annotate('{:.1f}'.format(sval), (r1,r2),  
                fontsize=9, ha='center', va='center')

            if smin:
                dotted_line = plt.Line2D((x1, x2), (y1, y2), lw=1., 
                                         ls='-', marker='.', 
                                         markersize=10, 
                                         markerfacecolor='r', 
                                         markeredgecolor='r', 
                                         alpha=0.7)

            else: 
                dotted_line = plt.Line2D((x1, x2), (y1, y2), lw=1., 
                                         ls='-', marker='.', 
                                         markersize=10, 
                                         markerfacecolor='r', 
                                         markeredgecolor='r', 
                                         alpha=0.3)

            plt.gca().add_line(dotted_line)
            """
            plt.annotate(r, (cx, cy), color='w', weight='bold', 
                fontsize=6, ha='center', va='center')
            """

        ## the two smallest
        iq= np.array(qvals).argsort()[:2]
        for k in range(len(qcontacts)):
            if k in iq: 
                connect_segments(qcontacts[k], qvals[k], True)
            else:
                connect_segments(qcontacts[k], qvals[k], False)


        colors = 100*np.random.rand(len(patches))
        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
        p.set_array(np.array(colors))
        ax.add_collection(p)


        pmax =  max(list(itertools.chain.from_iterable(peptides_1)))
        qmax =  max(list(itertools.chain.from_iterable(peptides_2)))


        M = mhcMod.MHCIModel()
        M.domain_architecture(fig, ax)

        T = tcrMod.TCRModel()
        T.domain_architecture(fig, ax)

        xmax = max([pmax,qmax])
        xmax=200
        
        ax.set_xlim([-10, xmax+10])
        ax.set_ylim([0, 10.5])
        



# ------------------------------------------
if __name__ == '__main__':

#'./res_3VXS/pairList.pkl
#'./res_3VXU/pairList.pkl'
#'./res_3KPR/pairList.pkl
    p=['./res_1LP9/pairList.pkl','./res_1MI5/pairList.pkl','./res_1OGA/pairList.pkl','./res_1QRN/pairList.pkl','./res_2BNQ/pairList.pkl','./res_2BNR/pairList.pkl','./res_2F53/pairList.pkl','./res_2F54/pairList.pkl','./res_2J8U/pairList.pkl','./res_2JCC/pairList.pkl','./res_2P5E/pairList.pkl','./res_2P5W/pairList.pkl','./res_2PYE/pairList.pkl','./res_2UWE/pairList.pkl','./res_2VLJ/pairList.pkl','./res_2YPL/pairList.pkl','./res_3D39/pairList.pkl','./res_3D3V/pairList.pkl','./res_3GSN/pairList.pkl','./res_3H9S/pairList.pkl','./res_3HG1/pairList.pkl','./res_3KPS/pairList.pkl','./res_3UTS/pairList.pkl','./res_3VXM/pairList.pkl','./res_3VXR/pairList.pkl','./res_3W0W/pairList.pkl','./res_4FTV/pairList.pkl','./res_4JFE/pairList.pkl','./res_4JFF/pairList.pkl','./res_4JRX/pairList.pkl','./res_4JRY/pairList.pkl','./res_4L3E/pairList.pkl','./res_4MNQ/pairList.pkl','./res_4PRI/pairList.pkl','./res_4PRP/pairList.pkl','./res_4QOK/pairList.pkl','./res_4QRR/pairList.pkl']


    #p = ['./res_4FTV/pairList.pkl']

    #clade 2
    p=['./res_1LP9/pairList.pkl','./res_2J8U/pairList.pkl', './res_2JCC/pairList.pkl','./res_2UWE/pairList.pkl']

    #clade 4
    p=['./res_4PRP/pairList.pkl']

    #clade 7
    p=['./res_3VXM/pairList.pkl']

    #clade 11
    p=['./res_3KPS/pairList.pkl', './res_1MI5/pairList.pkl']

    #clade 13

    #clade 14
    p=['./res_4JRX/pairList.pkl']



    #clade 16  
    p =[] 

    #clade 17
    p = [ './res_2BNR/pairList.pkl', './res_2PYE/pairList.pkl', './res_2BNQ/pairList.pkl', './res_2P5W/pairList.pkl', './res_2F53/pairList.pkl','./res_2F54/pairList.pkl', './res_3VXR/pairList.pkl']


    #clade 21
    p = ['./res_4FTV/pairList.pkl', './res_3D39/pairList.pkl', './res_3D3V/pairList.pkl', './res_3UTS/pairList.pkl','./res_3HG1/pairList.pkl']

    for kp in p: 
        #C = trMH.mhcComplex( 'res_4FTV/pairList.pkl' )
        print "--------------------"
        print kp
        C = trMH.mhcComplex( kp )
        L = C.MhcComplexList
        for rho in range(len(L)):
            print L[rho].pairs
            bStruct =  L[rho]

            D = plotPeptideMap ( kp.split("/")[1], bStruct )



