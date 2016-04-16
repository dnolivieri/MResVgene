#!/usr/bin/env python

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Motif
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp


query="GLSQVQLVQSEAELKKPGISVKVYCKASGYIFTNYNMHWVRQIPGHGFEWIGVISPSSGGTSYAQKFQGRVTMTRETSTSTAYMDLSSLGSEDIAMY"

query="RTREGGPGNAIEGQVRCAIRGRADNQGGLDVLSGQGPLWGLPWLWLVGVGALSGAAVLVVGGVLGALSAAVLSVAEGGLGVLSEAIVLADKLG"


blastp_cmdline = Blastp(db="./vsCdb/vsCdb", outfmt=6)

stdout, stderr = blastp_cmdline(stdin=query)

print "stdout=", stdout
print "stderr=", stderr
#s = stdout.split("\t")[10]

if stdout==None:
    print "fuck"

if stdout=="":
    print "off"


"""
bl_score = float(s)
print bl_score, type(bl_score)


if bl_score < 1.0e-3: 
    print "success"
"""
