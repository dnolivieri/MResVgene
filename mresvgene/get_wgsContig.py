#!/usr/bin/env python
"""
   dnolivieri:  updated ...4 nov 2015
      *** Just gets the contigs specified....

"""

import numpy as np
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
import re
import json
import cPickle as pickle

import errno

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


class GetWGSContig:
    def __init__(self, S):
        self.S = S

    def get_wgs_contigs(self, pList, infile, sdir):
        ofile_name = sdir +"/"+ os.path.basename(infile)
        ofile = open(ofile_name, "w")
        rec_cnt=0
        for record in SeqIO.parse(infile, "fasta"):
            if record.name.split("|")[1]  in pList:
                print record.name
                SeqIO.write(record ,ofile, "fasta")
                rec_cnt+=1

        ofile.close()

    def run(self):
        sbar= self.S["pogona"]["wgs"]
        sdir= self.S["pogona"]["sdir"]
        print sbar
        print sdir
        pList = ['CEMB01007043','CEMB01006509']

        self.get_wgs_contigs(pList, sbar, sdir)


# --------------------------------------
if __name__ == '__main__':

    WGS_files = 'WGS_files.json'
    json_data=open(WGS_files)
    S = json.load(json_data)
    json_data.close()

    
    R = GetWGSContig( S )
    R.run()



