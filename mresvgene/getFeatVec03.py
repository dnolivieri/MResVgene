#!/usr/bin/env python
"""
dnolivieri: (started: 9 sept 2014)
   modified for the MHC.
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
from Bio import SeqIO
from Bio.Seq import Seq

from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import json
import cPickle as pickle
import timeit

rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

AA = {0:'A',1:'R',2:'N',3:'D',4:'C',5:'Q',6:'E',7:'G',8:'H',9:'I',10:'L',11:'K',12:'M',13:'F',14:'P',15:'S',16:'T',17:'W',18:'Y',19:'V'}



class GetFeatVectors:
    def __init__(self, S,  desc_method):
        self.S = S
        self.desc_method   = desc_method
        qp = open('./mats/normalizedAA_Matrix.pkl', 'rb')
        self.normI = pickle.load(qp)

        #self.analyze_files()
        #self.get_background()


    def getPDT3(self, seq):
        Dvec={}
        cnt=0
        for q in self.normI:
            sumDseq=0.0
            for i in range(len(seq)-1):
                pi1 = rno[seq[i]]
                pi2 = rno[seq[i+1]]
                sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
                #sumDseq+= (q[rno[seq[i]]] - q[rno[seq[i+1]]])*(q[rno[seq[i]]] - q[rno[seq[i+1]]])

            sumDseq = sumDseq/np.float(len(seq)-1)
            Dvec.update( {str(cnt): sumDseq} )
            cnt+=1

        return Dvec

    
    def getPDT3B(self, seq):
        Dvec={}
        cnt=0
        for q in self.normI:
            for qLamda in range(1,4):
                sumDseq=0.0
                for i in range(len(seq)- qLamda):
                    pi1 = rno[seq[i]]
                    pi2 = rno[seq[i+qLamda]]
                    sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
                sumDseq = sumDseq/np.float(len(seq)-1)
                Dvec.update( {str(cnt): sumDseq} )
                cnt+=1
        return Dvec



    def get_background(self, Nmax, iterCount ):
        D=self.total_random_backgrnd(Nmax)
        npzoutfile="./bstrap/train_bckgndSignal_"+str(iterCount)+ ".npz"
        np.savez(npzoutfile, dp=np.array(D) )

        return np.array(D)


    def total_random_backgrnd(self, Nmax):
        qbar=[]
        for i in range(Nmax):
            indx = list(np.random.randint(0,19, 92))
            seqList = [ AA[j] for j in indx ]
            #print seqList
            seq =''
            rbar = seq.join(seqList)
            descObject=PyPro.GetProDes(rbar)
            if self.desc_method=='AAComp':
                T = descObject.GetAAComp()
            elif self.desc_method=='GestQSO':
                T = descObject.GetQSO()
            elif self.desc_method=='GetGearyAuto':
                T = descObject.GetMoreauBrotoAuto()
            elif self.desc_method=='GetCTD':
                T=descObject.GetCTD()
            elif self.desc_method=='GetPAAC':
                T =descObject.GetPAAC()
            elif self.desc_method=='PDT':
                T = self.getPDT3(rbar)
            elif self.desc_method=='PDT2':
                print "PDT3B"
                T = self.getPDT3B(rbar)
            else:
                print "here"
                T=descObject.GetCTD()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print i, Tx[0:5]
            #raw_input('press to continue')
            qbar.append(Tx)

        return np.array(qbar)



    def exon1_random_bkg(self, Nmax):
        qbar=[]
        cnt=0
        for i in range(Nmax):
            indx = list(np.random.randint(0,19, 92))
            seqList = [ AA[j] for j in indx ]
            #print seqList
            while ('S' in seqList[0:3]) or ('H' in seqList[0:3]) or (not 'C' in seqList[5:20]) or ('G' in seqList[19:27]) or ('D' in seqList[26:30]) or (not 'G' in seqList[1:5]) or (not 'R' in seqList[4:6]) or (not 'E' in seqList[0:4]): 
                indx = list(np.random.randint(0,19, 92))
                seqList = [ AA[j] for j in indx ]
                #print "......seqList=", seqList

            """
            lx=np.random.randint(0,11)            
            print "lx=", lx
            if lx<2:
                x = np.random.randint(0,10)
                if not 'C' in seqList[x:x+2]:
                    seqList[x:x+2] = ['G','G']
                print "seqList=", seqList
                kx=np.random.randint(0,11)
                if kx<4:
                    x = np.random.randint(20,50)
                    seqList[x:x+2] = ['K','K']

            kx=np.random.randint(0,11)            
            print "kx=", kx
            if lx>9 and kx<4:
                x = np.random.randint(30,50)
                seqList[x:x+2] = ['K','K']

            kx=np.random.randint(0,11)            
            print "kx=", kx
            if kx<1:
                x = np.random.randint(30,50)
                seqList[x:x+2] = ['E','E']
            
            kx=np.random.randint(0,11)            
            print "kx=", kx
            if kx<1:
                x = np.random.randint(40,60)
                seqList[x:x+2] = ['S','S']

            """  

            #print seqList
            seq =''
            rbar = seq.join(seqList)
            descObject=PyPro.GetProDes(rbar)
            if self.desc_method=='AAComp':
                T = descObject.GetAAComp()
            elif self.desc_method=='GestQSO':
                T = descObject.GetQSO()
            elif self.desc_method=='GetGearyAuto':
                T = descObject.GetMoreauBrotoAuto()
            elif self.desc_method=='GetCTD':
                T=descObject.GetCTD()
            elif self.desc_method=='GetPAAC':
                T =descObject.GetPAAC()
            elif self.desc_method=='PDT':
                T = self.getPDT3(rbar)
            elif self.desc_method=='PDT2':
                T = self.getPDT3B(rbar)
            else:
                print "here"
                T=descObject.GetCTD()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print i, Tx[0:5]
            #raw_input('press to continue')
            qbar.append(Tx)
        return np.array(qbar)


    def exon2_random_bkg(self, Nmax):
        qbar=[]
        cnt=0
        for i in range(Nmax):
            indx = list(np.random.randint(0,19, 92))
            seqList = [ AA[j] for j in indx ]
            while ('H' in seqList[1:4]) or  ('GC' in seqList[5:27]) or  ('D' in seqList[25:32]) or ('C' in seqList[70:82]):
                indx = list(np.random.randint(0,19, 92))
                seqList = [ AA[j] for j in indx ]
                #print "......seqList=", seqList

            seq =''
            rbar = seq.join(seqList)
            descObject=PyPro.GetProDes(rbar)
            if self.desc_method=='AAComp':
                T = descObject.GetAAComp()
            elif self.desc_method=='GestQSO':
                T = descObject.GetQSO()
            elif self.desc_method=='GetGearyAuto':
                T = descObject.GetMoreauBrotoAuto()
            elif self.desc_method=='GetCTD':
                T=descObject.GetCTD()
            elif self.desc_method=='GetPAAC':
                T =descObject.GetPAAC()
            elif self.desc_method=='PDT':
                T = self.getPDT3(rbar)
            elif self.desc_method=='PDT2':
                T = self.getPDT3B(rbar)
            else:
                print "here"
                T=descObject.GetCTD()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print i, Tx[0:5]
            #raw_input('press to continue')
            qbar.append(Tx)

        return np.array(qbar)


    def exon3_random_bkg(self, Nmax):
        qbar=[]
        cnt=0
        for i in range(Nmax):
            indx = list(np.random.randint(0,19, 92))
            seqList = [ AA[j] for j in indx ]
            while ('C' in seqList[17:23]) or ('L' in seqList[13:20]) or ('F' in seqList[19:25]) or ('P' in seqList[22:28]) or ('C' in seqList[70:92]) or ('W' in seqList[80:92]):
                indx = list(np.random.randint(0,19, 92))
                seqList = [ AA[j] for j in indx ]


            seq =''
            rbar = seq.join(seqList)
            descObject=PyPro.GetProDes(rbar)
            if self.desc_method=='AAComp':
                T = descObject.GetAAComp()
            elif self.desc_method=='GestQSO':
                T = descObject.GetQSO()
            elif self.desc_method=='GetGearyAuto':
                T = descObject.GetMoreauBrotoAuto()
            elif self.desc_method=='GetCTD':
                T=descObject.GetCTD()
            elif self.desc_method=='GetPAAC':
                T =descObject.GetPAAC()
            elif self.desc_method=='PDT':
                T = self.getPDT3(rbar)
            elif self.desc_method=='PDT2':
                T = self.getPDT3B(rbar)
            else:
                print "here"
                T=descObject.GetCTD()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print i, Tx[0:5]
            #raw_input('press to continue')
            qbar.append(Tx)

        return np.array(qbar)




    def analyze_files(self):
        for loci in self.S.iterkeys():
            print loci
            test_file= self.S[loci]
            D = self.descriptors_from_fasta(test_file)
            #print "len(D)=", len(D)

            npzoutfile=test_file.replace("fasta","npz")
            np.savez(npzoutfile, dp=np.array(D) )

    def descriptors_from_fasta(self, infile):
        sbar=[]
        qbar=[]
        cnt=0
        for record in SeqIO.parse(infile, "fasta"):
            descObject=PyPro.GetProDes(record.seq.tostring())
            if ('X' not in record.seq.tostring()) and ('Z' not in record.seq.tostring()) and ('B' not in record.seq.tostring()):
                if self.desc_method=='AAComp':
                    T = descObject.GetAAComp()
                elif self.desc_method=='GestQSO':
                    T = descObject.GetQSO()
                elif self.desc_method=='GetGearyAuto':
                    T = descObject.GetMoranAuto()
                elif self.desc_method=='GetCTD':
                    T=descObject.GetCTD()
                elif self.desc_method=='GetPAAC':
                    T =descObject.GetPAAC()
                elif self.desc_method=='PDT':
                    T = self.getPDT3(record.seq.tostring())
                elif self.desc_method=='PDT2':
                    T = self.getPDT3B(record.seq.tostring())
                else:
                    T=descObject.GetCTD()
                Tx = [ T[x]  for x in T.iterkeys() ]
                #print Tx
                #raw_input('press to continue')
                print cnt, Tx[0:5]
                qbar.append(Tx)
                cnt+=1
                if cnt>1e9:
                    break
        return np.array(qbar)



## ---------------MAIN ----------------------------------
if __name__ == '__main__':


    Vs_testing = 'MHC_test.json'
    json_data=open( Vs_testing )
    S = json.load(json_data)
    json_data.close()

    B =  GetFeatVectors(S,  desc_method='PDT2')
    B.analyze_files()
