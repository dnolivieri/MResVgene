#!/usr/bin/env python
"""
### This is to get info: 
grep \> all_R.fasta | awk -F"|" '{split($1,a,"-"); s=sprintf("%s-%s", substr(a[2],1,6),a[3]); locus[s]+=1}END{for (i in locus){split(i,b,"-"); print b[2],i, locus[i]}}' | sort

"""
import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import BarPlot01 as Bplt

class LociRSdata:
    def __init__(self):

        self.taxa = {
            "AQIB01": "Chlorocebus",
            "CABD02": "Gorilla_gorilla", 
            "AQIA01": "Macaca_fascicularis",
            "JYKQ01": "Mandrillus_leucophaeus",
            "ABDC01": "Microcebus_murinus",
            "ADFV01": "Nomascus_leucogenys",
            "AJFE01": "Pan_paniscus",
            "AHZZ01": "Papio_anubis",
            "ABGA01": "Pongo_abelii",
            "JZKE01": "Propithecus_coquereli",
            "JABR01": "Rhinopithecus_roxellana",
            "AGCE01": "Saimiri",
            "ABRT01": "Tarsius_syrichta"
            }
            
    def get_data(self, fR, fS): 
        def parse_file (fbar):
            rfile= open(fbar,"r")
            d={}
            for k in self.taxa.keys():
                d.update( {k:{'ighv':0, 'iglv':0,'igkv':0, 'trav':0, 'trbv':0, 'trdv':0, 'trgv':0} }  )


            print d
            for l in rfile:
                lb=l.split()
                if lb!=[]: 
                    species=lb[1].split("-")[0]
                    locus=lb[0]
                    nseqs=lb[2]
                    print species, locus, nseqs
                    d[species].update({locus:nseqs})


            rfile.close()
            return d

        dR = parse_file(fR)
        dS = parse_file(fS)
        return dR, dS



    def make_plots(self,fR, fS):
        dR, dS = self.get_data(fR, fS)
        print "dR=", dR
        print "-----------------"
        print "dS=", dS
        return dR, dS

    def run(self):
        species=["AQIB01","CABD02","AQIA01","JYKQ01","ABDC01",
                 "ADFV01","AJFE01","AHZZ01","ABGA01","JZKE01",
                 "JABR01","AGCE01", "ABRT01"]
        loci= ['ighv','iglv', 'igkv','trav','trbv','trdv','trgv']

        dR,dS = self.make_plots(fR, fS)        
        M = Bplt.MakeBarPlots(dR, dS, species, loci)
        V,W = M.get_barvectors()
        M.barplot(V, W)

# -----------------------------------------------
if __name__ == '__main__':

    #fR = "./trees_orig/res_all_R.dat"
    #fS = "./trees_orig/res_all_S.dat"

    fR = "./trees/res_all_R.dat"
    fS = "./trees/res_all_S.dat"

    L = LociRSdata()
    L.run()


