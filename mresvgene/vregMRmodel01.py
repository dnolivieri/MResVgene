#!/usr/bin/env python
"""
   dnolivieri:  updated ...8 december 2015
      vregMRmodel

      The specific multi-resolution (MR) model for V-regions.


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
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json
import cPickle as pickle
from collections import defaultdict
from banyan import *
import multiprocessing
from copy import deepcopy

import timeit
import operator



def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)


class VregMRmodel:
    def __init__(self, iterCount, loci_classes, adapt_threshold):
        self.iterCount = iterCount
        self.adapt_threshold = adapt_threshold
        self.loci_classes= loci_classes
        self.rec_id=None
        self.rec_name=None
        self.rec_description=None
        self.outfile = None
        self.exfile=None
        self.nlevel = 2
        self.rfmodels = self.get_models()

       
    def set_record(self, rec_id, rec_name, rec_description):
        self.rec_name=rec_name
        self.rec_id=rec_id
        self.rec_description=str(rec_description).replace(",", " ")
        print "in set_record:", rec_name

       
    def set_outfile(self, outfile):
        self.outfile=outfile

    def set_exon_outfiles( self, exfile):
        self.exfile=exfile

        
    def get_models(self):
        rfmodels = []
        for loci in self.loci_classes: 
            nl=[]
            for mrlevel in range(np.power(2,self.nlevel)-1):
                matfile = "./bstrap/trainMat_" + loci + "_r"+str(mrlevel)+"_n"+ str(self.iterCount)+".pkl"
                print matfile
                fp = open(matfile, 'rb')
                nl.append( pickle.load(fp) )

            rfmodels.append( nl )
        #print rfmodels
        return rfmodels

    def maxLk_interval(self, z, zs): 

        print "****----in maxLk -----"
        print "z=", z,  ".....argmax=", np.argmax(z)
        print "zs=", zs, ".....argmin=", np.argmin(zs)
        print "***-....."

        izmax = np.argmax(z)
        zmax  = np.max(z)

        print "izmax=", izmax
        print "zmax=", zmax



        rb = [k for k in range(len(z)) if np.abs(zmax - z[k]) < 0.025  ]
        if len(rb)>1: 
            sb = [ zs[k] for k in rb ]
            izmax = np.argmin(sb)


        print "final izmax=", izmax
        return izmax


    def exon_MRprobabilities(self, reduced_seqbuffer):
        qbar =  self.sequence_MRhomology(reduced_seqbuffer)
        sbar=[]
        for q in qbar:
            ## q[0] positions; q[2] kmax; q[3] max prob;  q[4] prob MRlevels; q[5] sum(probMR_levels)
            s=(q[0],q[2],q[3], q[4], q[5] )
            sbar.append(s)            
            print "q=",q
            
        ExonList = self.descriminate_by_exon(sbar)           
        return ExonList


    def sequence_MRhomology(self, reduced_seqbuffer):
        # returns the sequences having homology at MR level.
        pos_seqindx = []
        LprobX=[]
        seqs=[]
        qseq = ()

        ## Step 1:  Get data from reduced sequence buffer:
        seqpos =  [x[0] for x in reduced_seqbuffer if len(x[2])!=0 ]
        seqbuffer =  [ str(x[1]) for x in  reduced_seqbuffer if len(x[2])!=0   ]
        Y  =  np.array([ x[2] for x in reduced_seqbuffer if len(x[2])!=0] )

        exonType=''
        print "Y=", Y.shape

        ## Step 2:  Reform matrix, Produce a numpy array for each MR level
        Ybar= []
        for k in range(np.power(2,self.nlevel)-1):
            Ybar.append(np.array([list(itertools.chain(*x[k])) for x in Y ]))


        ## Step 3: For each MultRes level matrix, run prediction with approp model 
        for loci in range(len(self.loci_classes)): 
            probX=[]
            for mrlevel in range(np.power(2,self.nlevel)-1):
                probX.append( self.rfmodels[loci][mrlevel].predict_proba(Ybar[:][mrlevel]) )
                #print type(Ybar[mrlevel]), Ybar[mrlevel].shape, Ybar[mrlevel][0:5]

            LprobX.append(probX)

        #print LprobX

        ## Step 4: Obtain Dict that summarizes data from each loci and all candidates 
        D={}
        for loci in self.loci_classes:
            D.update({loci:[]})
        for lk in range(len(self.loci_classes)): 
            loci = self.loci_classes[lk]
            probX= LprobX[lk]
            for j in range(probX[0].shape[0]):
                # since proba is symmetric, only keep positive val...
                Xbar = np.array([ probX[k][j][1] for k in range(len(probX)) ])
                print "Xbar=", Xbar
                D[loci].append(Xbar)
        print D

        ## Step 5: The logic for selecting Maximums:
        """
           - at the moment use a simple logic; just 
             sum the positive class probabilities from each MR.
             In this way, the one with the highest value wins!..
        """

        print "STEP 5 ---------------------------"
        for j in range( len(D[self.loci_classes[0]])):  
            rho = [ D[x][j] for x in D.iterkeys() ]
            tau = np.array([ D[x][j].sum() for x in D.iterkeys() ])
            loci_kmax = np.argmax(tau)
            loci_ymax = np.max(tau)
            print "rho=", rho, "    tau=", tau, ".....", loci_kmax, loci_ymax

            if ( np.max(rho[loci_kmax])  > 0.4 ):
                print "found ....", rho[loci_kmax]
                Pr_outclass = rho[loci_kmax]
                seqs.append( ( seqpos[j],seqbuffer[j], loci_kmax, np.max(rho[loci_kmax]), Pr_outclass, np.sum(Pr_outclass) ) )


        return seqs


    def descriminate_by_exon(self, Q):
        print "Q=", Q
        ## Q[0] positions; Q[1] kmax; Q[2] max prob;  q[3] Pr_outclass
        s = [x[0] for x in Q]
        eprm = [x[1] for x in Q]
        score= [x[2] for x in Q]
        pclass = [ x[3] for x in Q] 
        sclass = [ x[4] for x in Q] 


        print "s=", s
        print "eprm=", eprm
        print "score=", score
        print "pclass=", pclass
        print "sclass=", sclass
        
        tree =  SortedSet( s, key=lambda (start, end): (end, (end - start)),updator=OverlappingIntervalsUpdator)

        print tree
        r = set([])
        tbar=tree.copy()
        for sprm in s:
            a=tree.overlap(sprm)
            b=[s.index(x) for x in tree.overlap(sprm) ]
            c=[eprm[x] for x in b ]
            d=[score[x] for x in b ]

            pbar = [pclass[x] for x in b]
            spbar = [sclass[x] for x in b]

            gX=sorted(list_duplicates(c))
            if len(gX)>0:
                g=sorted(list_duplicates(c))[0][1]
                z=[d[i] for i in g ]
                z=   [d[i] for i in g ]
                zp = [pbar[i] for i in g]
                zs = [spbar[i] for i in g]

                #print "z=", z,  ".....argmax=", np.argmax(z)
                #print "zp=", zp, ".... zs=", zs, ".....argmin=", np.argmin(zs)
                izmax = self.maxLk_interval( z, zs )

                for i in g:
                    #if i!= g[np.argmax(z)]:
                    if i!= g[izmax]: 
                        tree.remove( s[b[i]] )
                        r.update([b[i]])
                else:
                    z=[]

        r_list = [i for j, i in enumerate(Q) if j not in list(r)]
        c=[ x[1] for x in r_list ]

        # ---- IN CLASS Overlaps ----------------

        tmpbuff=[]
        rindx_buffer=[]
        for i in range(1,len(c)):
            if (i-1) not in tmpbuff:
                tmpbuff.append(i-1)
            if c[i]== c[i-1]:
                if i not in tmpbuff:
                    tmpbuff.append(i)
            else:
                if len(tmpbuff)>0:
                    z = [ r_list[x][2] for x in tmpbuff ]
                    imax=np.argmax(z)
                    jmax=tmpbuff[imax]
                    for j in tmpbuff:
                        if j!=jmax:
                            rindx_buffer.append(j)

                tmpbuff=[]
            lastI=i

        merged = rindx_buffer
        Q = [i for j, i in enumerate(r_list) if j not in merged]

        ## ---- Out-Of-Class Overlaps -------------------------

        s = [ x[0] for x in Q ]
        eprm = [x[1] for x in Q]
        score= [x[2] for x in Q]

        tree =  SortedSet( s, key=lambda (start, end): (end, (end - start)),updator=OverlappingIntervalsUpdator)
        r = set([])
        for sprm in s:
            a=tree.overlap(sprm)
            b=[s.index(x) for x in tree.overlap(sprm) ]
            c=[eprm[x] for x in b ]
            d=[score[x] for x in b ]

            if len(b)>1:
                for i in range(len(d)):
                    if i!= np.argmax(d):
                        tree.remove( s[b[i]] )
                        r.update([b[i]])
                    else:
                        z=[]

        cand_list = [i for j, i in enumerate(Q) if j not in list(r)]
        print "cand_list=", cand_list
        return cand_list


    def V_exon_model( self, mcnt, seq, strand, AllExons):
        print "len(AllExons)=", len(AllExons)
        print "AllExons=", AllExons
        VList=[]


        ExonList = []
        for p in AllExons:
            indx1 = p[0]
            loci_index = p[1]
            rfscore1= p[2]
            exon1 = seq[indx1[0]+2: indx1[0]+ (indx1[1] - indx1[0])]
            #print   "----exon1--------------", "len(exon1)=", len(exon1)
            #print "start=",indx1[0]+2, "stop=", indx1[0]+ (indx1[1] - indx1[0])
            #print "exon1=",exon1
            #print "***exon1AA=",exon1[2:].translate()
            #print  seq[indx1[0]: indx1[0]+ 10], ".....", seq[ indx1[0]+(indx1[1] - indx1[0])-10:  indx1[0]+(indx1[1] - indx1[0])+2 ]
            #print  "  ",seq[indx1[0]+2: indx1[0]+ 10], ".....", seq[ indx1[0]+(indx1[1] - indx1[0])-10:  indx1[0]+(indx1[1] - indx1[0]) ]
            RNA = exon1
            Vreg= RNA[2:-1].translate()

            print "Vreg=", Vreg
            exon1_ival = (indx1[0], indx1[1])
            exons_rfscore =  rfscore1 

            ExonList.append( (RNA, self.rec_id, self.rec_description, exon1_ival, rfscore1, loci_index) )
            VList.append( (Vreg, self.rec_id, self.rec_description, exon1_ival, exons_rfscore, loci_index) )


        tmpcnt = mcnt
        excnt=tmpcnt
        for p in ExonList:
            print "exon1=", str(p[0])
            rec_id = p[1].split("|")[2]
            t= p[2].split("|")
            rec_name=t[0]+"_"+t[1]
            locus= self.loci_classes[p[5]]
            recordB=SeqRecord(p[0], id = "exon1-"+str(excnt)+"-RF-"+rec_id+"|"+rec_name+"|"+locus, description=str(p[3])+"|"+str(p[4])+"|"+ str(strand))
            self.exfile.write(recordB.format("fasta"))
            excnt+=1

        for p in VList:
            print "OUTPUT V=", str(p[0])
            rec_id = p[1].split("|")[2]
            t= p[2].split("|")
            rec_name=t[0]+"_"+t[1]
            locus= self.loci_classes[p[5]]
            recordB=SeqRecord(p[0], id = "V-"+str(mcnt)+"-RF-"+rec_id+"|"+rec_name+"|"+locus, description=str(p[3])+"|"+str(p[4])+"|"+ str(strand))
            self.outfile.write(recordB.format("fasta"))

            mcnt+=1


        return mcnt
        
#-------------------
if __name__ == '__main__':

    iterCount=0
    adapt_threshold = 0.6    
    loci_classes= ['ighv']

    V = VregMRmodel(iterCount, loci_classes, adapt_threshold)
