#!/usr/bin/env python
"""
   dnolivieri:  updated ...26 december 2015
      - modified to use the boostrap code.

      This does the prediction of Vs based 
      upon a multi-resolution training.

      - also should implement (eventually) the 
       bootstrap learning.
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
import vregMRmodel04 as VRmodel
from MResStruct02 import MResStruct
from PDT3method import PDT3

# Note this redefinition to handle overlaps
#import regex as re
import re



class SimpleMapReduce(object):
    def __init__(self, map_func, partition_func, reduce_func, num_workers=None):
        self.map_func = map_func
        self.partition_func = partition_func
        self.reduce_func = reduce_func
        #self.pool = multiprocessing.Pool(num_workers)
        self.pool = multiprocessing.Pool(4)

    def partition(self, mapped_values):
        partitioned_data= self.partition_func( list(mapped_values))
        return partitioned_data

    def __call__(self, inputs, chunksize=1):
        map_responses = self.pool.map(self.map_func,inputs, chunksize=chunksize)
        partitioned_data = self.partition(itertools.chain(*map_responses))
        reduced_values = self.pool.map(self.reduce_func, partitioned_data)
        return reduced_values


# --------------HERE  IS THE Implementation ------

def divide_work(seq):
    def NNchunks(l, n, delta):
        for i in xrange(0, len(l), n):
            if i==0:
                yield (i, i+n,  l[i:i+n])
            else:
                yield (i-1000, i+n,  l[i-1000:i+n])

    return NNchunks(seq, 22000, 1000)




def process_blockNN(d):
    #print multiprocessing.current_process().name , d[0], d[1]

    sbar=d[2]
    ix=d[0]
    """
    ###if used with regex
    x=[i.start()+ix for i in re.finditer(r'AG', str(sbar), overlapped=True)]
    y=[i.start()+ix for i in re.finditer(r'CAC', str(sbar), overlapped=True)]
    """
    ### if used with re
    x=[i.start()+ix for i in re.finditer("AG", str(sbar))]
    y=[i.start()+ix for i in re.finditer("CAC", str(sbar))]


    s=[(i,j) for i,j in itertools.product(x,y) if j>i and ( np.abs(i-j)>270  and np.abs(i-j)<330) ]
    #print "---before translate"
    #print "y=", y
    #print "s=",s


    cnt=0
    cand_list=[]
    cand_seqs=[]
    for i in range(len(s)):
        idx1=s[i][0]-ix   # corrects block start
        idx2=s[i][1]-ix
        #test = sbar[idx1+4:idx1+(idx2-idx1)-1]   ## add -1
        test = sbar[idx1+4:idx1+(idx2-idx1)]   ## do not add zero
        #test = test[0:len(test)-(len(test)%3)]
        # (6jan2016): put some condition on seq
        has_delta=False

        ## this corrects the problems with Biopython inserting extra aa.
        delta = len(test)%3
        if delta!=0:
            test = test[:-delta]
            has_delta = True

        p = test.translate(to_stop=True)
        print i, s[i]
        print " test=",test, "....len(test)=", len(test)
        print " p=", p
        print " condition=", np.abs(int(len(test)/3.)-len(p))

        ## (10-march) changed condition to "<=1" because stop codon could occur at end
        if (np.abs(int(len(test)/3.)-len(p))<=20) and ( (len(p)>65) and (len(p) < 130)):
            if has_delta ==True: 
                s[i] = (s[i][0], s[i][1]-delta)

            cand_list.append(s[i])
            cand_seqs.append(p)
            print "-------------------"
            print s[i], ix, idx1, idx2, len(test)
            print "seq=", sbar[idx1:idx1+(idx2-idx1)+3]
            print "AA=", p
            print "-------------------" 
            print 
            cnt+=1

    output_list = zip(cand_list, cand_seqs)
    return output_list



def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)


def partition_func(s):
    #print "----start partition_func"
    #print "s=", s
    D = {}
    for x in s:
        D.update({x[0]:[]})

    for x in s:
        D[x[0]].append(x[1])

    cand_list = []
    for k in D.iterkeys():
        cand_list.append((k,D[k][0]))
    
    cand_list.sort()

    #print "---after, but from partition func---"
    #print "cand_list=",cand_list
    return cand_list


rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,
       'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,
       'J':10, 'B':2 }


qp = open('./aaRefs/normalizedAA_Matrix.pkl', 'rb')
normI = pickle.load(qp)
nlevel=2
Mstruct = MResStruct(nlevel)

def getPDT3(seq):
    Dvec={}
    cnt=0
    for q in normI:
        sumDseq=0.0
        for i in range(len(seq)-1):
            pi1 = rno[seq[i]]
            pi2 = rno[seq[i+1]]
            sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
        sumDseq = sumDseq/np.float(len(seq)-1)
        Dvec.update( {str(cnt): sumDseq} )
        cnt+=1
    return Dvec

def getPDT3B(seq):
    Dvec={}
    cnt=0
    for q in normI:
        for qLamda in range(1,4):
            sumDseq=0.0
            for i in range(len(seq)- qLamda):
                pi1 = rno[seq[i]]
                pi2 = rno[seq[i+qLamda]]
                sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
            sumDseq = sumDseq/np.float(len(seq)-1)
            Dvec.update( {str(cnt): sumDseq} )
            cnt+=1
    #print "len(Dvec)=",len(Dvec)
    return Dvec


def get_lambda_seq2vector(seq, qLambda):
    Dvec={}
    cnt=0
    for q in normI:
        sumDseq=0.0
        for i in range(len(seq)- qLambda):
            pi1 = rno[seq[i]]
            pi2 = rno[seq[i+qLambda]]
            sumDseq+= (q[ pi1 ] - q[ pi2 ])*(q[ pi1 ] - q[ pi2 ])
        sumDseq = sumDseq/np.float(len(seq)-1)
        Dvec.update( {str(cnt): sumDseq} )
        cnt+=1
    return Dvec



def count_sequential_doubles(seq):
    dbar = [ seq[i-1]+seq[i] for i in range(len(seq)) if i>0 ]
    res=collections.Counter(dbar)
    return res

def get_freq_seq2vector(seq):
    def merge_two_dicts(x, y):
        '''Given two dicts, merge them into a new dict as a shallow copy.'''
        z = x.copy()
        z.update(y)
        return z
    aa ='ARNDCQEGHILKMFPSTWYV'
    ## single
    dbar = [ x for x in aa] 
    S={}
    for i in dbar:
        S.update({i:0})
    T=collections.Counter(seq)
    y = merge_two_dicts(S, T)        
    ybar = np.array([ y[j] for j in dbar])

    ## double.
    dbar = [ x[0]+x[1] for x in list(itertools.product(aa, repeat=2))] 
    D={}
    for i in dbar:
        D.update({i:0})
    B= count_sequential_doubles(seq)
    z = merge_two_dicts(D, B)
    zbar = np.array([ z[j] for j in dbar ])
    
    rho = np.hstack( (ybar,zbar))
    Dvec={}
    for k in range(len(rho)):
        Dvec.update( {str(k): rho[k]} )
    
    return Dvec



# ......
def reduce_candidates(s):
    #print multiprocessing.current_process().name
    Tx=[]
    seq=str(s[1])

    """
    ## 9dec: this code is replaced by the code seq_to_MRdict( seq )
    if ('X' not in  seq) and ('Z' not in seq) and ('B' not in seq):
        T = getPDT3B(seq)  ## should be MR....
        Tx = [ T[x]  for x in T.iterkeys() ]
    Y = np.array(Tx)
    """

    Y={}
    #print "in reduce_candidates:   s=", s
    if ('X' not in  seq) and ('Z' not in seq) and ('B' not in seq):
        Y = seq_to_MRdict( seq ) 

    return (s[0], s[1], Y)


# ---------------- interface to the above.......
def seq_to_MRdict( seq_data ): 
    nlevel=2
    D = {}
    #qLambda=1  ### this can change!!
    for k in range(np.power(2,nlevel)-1):
        D.update( {k:[]} )

    S  = Mstruct.get_pyramid(seq_data)
    knt=0
    for k in range(len(S)):
        qLambda = nlevel - k
        for kseq in S[k]:
            if k == 0:
                T =  get_freq_seq2vector(kseq)
            else:
                T =  get_freq_seq2vector(kseq)
                #T =  get_lambda_seq2vector(kseq, qLambda)
            
            Tx = [ T[str(x)]  for x in range(len(T)) ]
            D[knt].append(Tx)
            knt+=1
    return D

#----------------------------------
    
class RunVregPredict:
    def __init__(self, S,  speciesList, check_contigs, contigs=None):
        self.S = S
        self.speciesList = speciesList
        self.mapper = SimpleMapReduce(process_blockNN, partition_func, reduce_candidates)
        self.outDir = "./bstrap/"

        qp = open('./aaRefs/normalizedAA_Matrix.pkl', 'rb')
        self.normI = pickle.load(qp)
        self.predicted_seqs = []
        self.check_contigs = check_contigs
        """
        if check_contigs: 
            self.contigs = self.get_contigs(speciesList[0])
        else: 
            self.contigs = []
        """

        if contigs==None: 
            self.contigs = self.get_contigs(speciesList[0])
        else:
            self.contigs= contigs

        print "self.contigs=", self.contigs


    def get_contigs(self, species):
        contigs=[]
        fp = open(self.S[species]["cntgs"], "r")
        for lk in fp:
            contigs.append(lk.strip())
        fp.close()
        #print contigs
        return contigs


    def analyze_files(self, iterCount, loci_classes, adapt_threshold):
        Rmodel = VRmodel.VregMRmodel(iterCount, loci_classes, adapt_threshold)
        print "len(Rmodel.rfmodels)=", len(Rmodel.rfmodels)

        ofile = open("bkg_out.dat","a+")
        Rmodel.set_bckgoutfile( ofile )
        
        for species in self.speciesList:
            fbar= self.S[species]["WGS"]
            print fbar
            outFile = self.outDir + os.path.basename(fbar).replace(".fasta", "_"+str(iterCount)+"_outRF.fasta")
            ofile = open(outFile,"w")
            Rmodel.set_outfile( ofile )

            fb = self.outDir + os.path.basename(fbar).replace(".fasta", "_"+str(iterCount)+"_exon.fasta")  
            exfile1 = open(fb,"w")            
            Rmodel.set_exon_outfiles( exfile1 ) 

            start_time = timeit.default_timer()
            gene_cnt=0


            for strand in [1, -1]:
                qbar=deepcopy(self.contigs)
                print "STRAND=", strand
                for record in SeqIO.parse(fbar, "fasta"):
                    if self.check_contigs: 
                        if ( record.id.split("|")[3] not in self.contigs):
                            continue
                    print "record.id=", record.id
                    print "cnts=",record.id.split("|")[3]
                    print "qbar=", qbar
                    if self.check_contigs: 
                        qbar.remove(record.id.split("|")[3])
                    if strand == 1:
                        seq=record.seq
                    else:
                        seq=record.seq.reverse_complement()

                    Rmodel.set_record(record.id, record.name, record.description)
                    seq_size=len(seq)

                    res= self.mapper( divide_work(seq) )

                    """
                    print "len(res)=", len(res)
                    for ix in range(2):
                        print res[ix][0], res[ix][1], type(res[ix][2])
                    """

                    Elist=Rmodel.exon_MRprobabilities(res)
                    gene_cnt = Rmodel.V_exon_model(gene_cnt, seq, strand, Elist)
                    #res=None
                    #Elist=None
                    if len(qbar)==0: 
                        break

            ofile.close()
            elapsed = timeit.default_timer() - start_time
            print "ELAPSED TIME =", elapsed



        
#----------------------------------
if __name__ == '__main__':

    #Vs_prediction = 'vs_predict.json'
    WGS_files = 'WGS_files.json'
    json_data=open( WGS_files )
    S = json.load(json_data)
    json_data.close()

    #spList = ["pogona"]
    spList = ["Pan_paniscus"]
    print S
    M =  RunVregPredict(S,  speciesList= spList )

    iterCount=0
    adapt_threshold=0.6
    loci_classes=['ighv', 'igkv']
    M.analyze_files(iterCount, loci_classes, adapt_threshold)
