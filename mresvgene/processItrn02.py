#!/usr/bin/env python
"""
   dnolivieri: (27dec2015):
    - this does the iteration selection step for next iteration.
    - adapted from cleanIteration
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


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename



class KDEProbDistribution: 
    def __init__(self, cntIter):
        self.infile="./bstrap/bstrp_iteration_"+str(cntIter)+".fasta"
        self.cntIter = cntIter
        self.p=self.get_prob_fasta()



class ProcessIteration: 
    def __init__(self, p_iterN, loci_classes):
        self.p_iterN = p_iterN
        self.n_iterN = p_iterN + 1 
        self.loci_classes = loci_classes




    def write_out_full_iteration(self, fileList):
        #1. Get all sequences from this iteration and save in one file:
        outFile = "./bstrap/bstrp_iteration_"+str(self.p_iterN)+".fasta"
        ofile = open(outFile, "w")
        posCount=0
        for vFile in fileList:
            for record in SeqIO.parse(vFile, "fasta"):
                SeqIO.write(record, ofile, "fasta")
                posCount+=1

        ofile.close()
        print "posCount=", posCount

    def write_out_full_loci_species(self, fileList):
        for vFile in fileList:
            for loci in self.loci_classes: 
                outFile = vFile.replace(".fasta","_"+loci+".fasta")
                ofile = open(outFile, "w")
                for record in SeqIO.parse(vFile, "fasta"):
                    if record.id.split("-")[2]==loci: 
                        SeqIO.write(record, ofile, "fasta")
                ofile.close()

    def write_out_full_loci_iteration(self, fileList):
        print "fileList=", fileList
        for loci in self.loci_classes: 
            outFile = "./bstrap/"+loci +"_"+ str(self.p_iterN + 1) +".fasta"
            ofile = open(outFile, "w")
            for vFile in fileList: 
                for record in SeqIO.parse(vFile, "fasta"):
                    if record.id.split("-")[2]==loci: 
                        SeqIO.write(record, ofile, "fasta")


            ofile.close()


    def get_new_sequences(self):
        fileList = list(find_files("./bstrap/", "*_"+str(self.p_iterN) +"_outRF.fasta"))
        ## 1. all sequences output:        
        self.write_out_full_iteration(fileList)

        #2 write each loci.
        self.write_out_full_loci_iteration(fileList)

        # 3. first write the reference file: 
        for loci in self.loci_classes: 
            refFile = "./bstrap/"+loci + "_init.fasta"
            outFile = "./bstrap/"+loci +"_"+ str(self.p_iterN + 1) +".fasta"
            ofile = open(outFile, "w")
            for record in SeqIO.parse(refFile, "fasta"):
                SeqIO.write(record, ofile, "fasta")

            for vFile in fileList: 
                inFile = vFile.replace(".fasta","_"+loci+".fasta")                
                print loci, "inFile=", inFile 
                for record in SeqIO.parse(inFile, "fasta"):
                    SeqIO.write(record, ofile, "fasta")
            ofile.close()




    def get_difference_seqs(self):
        fileList = list(find_files("./bstrap/", "*_"+str(self.p_iterN) +"_outRF.fasta"))
        def get_tuple_dictionary( inFile ):
            T = {}
            for record in SeqIO.parse(inFile, "fasta"):
                cntg=record.id.split("-")[1]
                ptuple=record.description.split("|")[1]
                print cntg, ptuple
                if cntg not in list(T.iterkeys()):
                    T.update({cntg: [ptuple]} )
                else:
                    T[cntg].append(ptuple)
            return T

        deltaFileList=[]
        for vFile in fileList:
            Tj = get_tuple_dictionary( vFile )
            wFile = vFile.replace(str(self.p_iterN) +"_outRF", str(self.p_iterN-1) +"_outRF")
            print wFile
            Tjm1 = get_tuple_dictionary( wFile )

            tjkeys = list(Tj.iterkeys())
            tjm1keys = list(Tjm1.iterkeys())

            print tjkeys
            print tjm1keys
            
            ## does intersection
            common_keys = list(set(tjkeys) & set(tjm1keys))
            ## does difference (but only in tj)
            uniquej_keys = list(set(tjkeys) - set(tjm1keys))
            print common_keys
            print uniquej_keys

            DeltaPj = []
            ## first add the uniquej pos: 
            for x in uniquej_keys: 
                DeltaPj.append( Tj[x] )

            ## now run over each common key and get what is new:
            for x in common_keys: 
                pj = list(set(Tj[x]) - set(Tjm1[x]))
                if pj != []: 
                    DeltaPj.append( pj )
                print x, pj
                ## remove common from T
                if pj==[]: 
                    Tj.pop(x)
                else: 
                    Tj[x] = pj

            print DeltaPj
            Dp = list(itertools.chain.from_iterable(DeltaPj))
            print Dp

            print "Tj=", Tj

            ## For the moment, just write out the delta file for each.
            outFile = vFile.replace("_outRF","_delta_outRF")
            deltaFileList.append(outFile)
            ofile = open(outFile, "w")
            
            for record in SeqIO.parse(vFile, "fasta"):
                cntg=record.id.split("-")[1]
                ptuple=record.description.split("|")[1]
                print cntg, ptuple
                if cntg in list(Tj.iterkeys()): 
                    if ptuple in Tj[cntg]: 
                        SeqIO.write(record, ofile, "fasta")                        
            ofile.close()

        return  deltaFileList


    def write_full_delta_files(self, deltaFileList):
        outFile = "./bstrap/bstrp_iteration_delta_"+str(self.p_iterN)+".fasta"
        ofile = open(outFile, "w")
        for dfile in deltaFileList: 
            for record in SeqIO.parse(dfile, "fasta"):
                SeqIO.write(record, ofile, "fasta")            
        ofile.close()


    def write_delta_loci_species(self, deltaFileList):
        for vFile in deltaFileList: 
            for loci in self.loci_classes: 
                outFile = vFile.replace(".fasta","_delta_"+loci+".fasta")
                ofile = open(outFile, "w")
                for record in SeqIO.parse(vFile, "fasta"):
                    if record.id.split("-")[2]==loci: 
                        SeqIO.write(record, ofile, "fasta")
                ofile.close()


    def write_delta_loci_iteration(self, deltaFileList):
        for loci in self.loci_classes: 
            outFile = "./bstrap/"+loci +"_delta_"+ str(self.p_iterN + 1) +".fasta"
            ofile = open(outFile, "w")
            for vFile in deltaFileList: 
                for record in SeqIO.parse(vFile, "fasta"):
                    if record.id.split("-")[2]==loci: 
                        print record.id
                        SeqIO.write(record, ofile, "fasta")
            ofile.close()


    def run(self):
        
        fileList = list(find_files("./bstrap/", "*_"+str(self.p_iterN) +"_outRF.fasta"))

        ## 1. all sequences output:        
        self.write_out_full_iteration(fileList)
        #2 write each loci.
        self.write_out_full_loci_species(fileList)
        self.write_out_full_loci_iteration(fileList)



        # 3. process all delta files
        if self.p_iterN>0:
            deltaFileList=self.get_difference_seqs()
            self.write_full_delta_files(deltaFileList)
            self.write_delta_loci_species(deltaFileList)
            self.write_delta_loci_iteration(deltaFileList)
        else: 
            self.write_delta_loci_iteration(fileList)
            


# ------------------------------------------
if __name__ == '__main__': 

    presentN=0
    loci_classes=['ighv','igkv']
    P= ProcessIteration(presentN, loci_classes)
    #P.get_new_sequences()
    #P.get_difference_seqs()
    P.run()

