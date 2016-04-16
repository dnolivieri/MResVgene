#!/usr/bin/env python
"""
   dnolivieri:  updated ...5-jan-2016
     --inverse string

    >>> 'hello world'[::-1]
    'dlrow olleh'
     This is extended slice syntax. It works by doing [begin:end:step] 
     - by leaving begin and end off and specifying a step of -1, 
       it reverses a string.

     """
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import re
import time
import os, fnmatch
import sys
import itertools


#-------------------
if __name__ == '__main__':


    #(78167, 78474) 303
    seq= "GTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTGACGACTACATGGAGTGGGTCCGCCAGGCTCCAGGAAAGGGGCTGGAGTGGGTTGGACAAATTAATCCTAATGGGGGTACCACATTCCTCATGGATTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTTTATCTGCAAATTAACAGCCTGAAAATAGAGGACACGGCCGTGTATTACTGTACTAGA"
    #AA= VQCEVQLVESGGGLVQPGGSLRLSCAASGFTFSDDYMEWVRQAPGKGLEWVGQINPNGGTTFLMDSVKGRFTISRDNAKNTLYLQINSLKIEDTAVYYCTR



    seq= "GTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTGACGACTACATGGAGTGGGTCCGCCAGGCTCCAGGAAAGGGGCTGGAGTGGGTTGGACAAATTAATCCTAATGGGGGTACCACATTCCTCATGGATTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTTTATCTGCAAATTAACAGCCTGAAAATAGAGGACACGGCCGTGTATTACTGTACTAGA"

    sbar="TTTCTAACCAGGGTGTCTCTGTGTTCGCAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTGACGACTACATGGAGTGGGTCCGCCAGGCTCCAGGAAAGGGGCTGGAGTGGGTTGGACAAATTAATCCTAATGGGGGTACCACATTCCTCATGGATTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTTTATCTGCAAATTAACAGCCTGAAAATAGAGGACACGGCCGTGTATTACTGTACTAGACACACAATGAGGGGAGGTCAGTGTGAGCCCAGACACAAA"

    """
    ix=0
    x=[i.start()+ix for i in re.finditer("AG", str(sbar))]
    y=[i.start()+ix for i in re.finditer("CAC", str(sbar))]
    s=[(i,j) for i,j in itertools.product(x,y) if j>i and ( np.abs(i-j)>300  and np.abs(i-j)<324) ]
    print "---before translate"
    print "y=", y
    print "s=",s
    """




    #print "seq=", seq
    #print
    #print "iseq=", seq[::-1]
    fbar="Macacaf_AQIA01066387.fasta"
    fbar="Macacaf_AQIA01017118.fasta"

    for record in SeqIO.parse(fbar, "fasta"):
        print "record.id=", record.id
        seq=record.seq.reverse_complement()
        print seq
