#!/usr/bin/env python
"""
   dnolivieri:  updated ...29 dec 2015
     used to generate the background.

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class GenBackground:
    def __init__(self,infile):
        self.inFile = infile


    def get_dist_data(self):

        #outFile = self.inFile.replace(".fasta", "_r6500.fasta")
        outFile = self.inFile.replace(".fasta", "_m.fasta")
        ofile = open(outFile ,"w")

        bknt=0
        for record in SeqIO.parse(self.inFile, "fasta"):
            p = record.seq
            z = p[33:45]
            z2 = p[18:30]
            z3 = p[-12:-1]
            """
            if ('C' in z2) and ('W' in z):
                print bknt, p
            """
                
            if ('C' not in z2) and ('W' not in z) and bknt<6500:
                SeqIO.write(record ,ofile, "fasta")
                bknt+=1

        print "bknt=", bknt
        ofile.close()


    def run(self):
        self.get_dist_data()


## ---------------MAIN ----------------------------------
if __name__ == '__main__':
    
    infile = "./bstrap/trdv_2.fasta"
    infile = "./bckgnd_run.fasta"

    infile = "./dataVgeneDB/All_bkg.fasta"

    T = GenBackground(infile)
    T.run()
    

