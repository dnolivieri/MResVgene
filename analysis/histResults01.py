#!/usr/bin/env python
"""
grep \> all_R.fasta | awk -F"-" '{split($5,a," "); print a[1]}'>  res_probR.dat

grep \> all_S.fasta | awk '{split($2,a,"|"); print a[1]}' > res_probS.dat


"""
import numpy as np
import scipy as sp
import os

import seaborn as sns
import matplotlib.pyplot as plt


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
            d=[]
            rfile= open(fbar,"r")
            for l in rfile:
                d.append(float(l))
            rfile.close()
            return d

        dR = parse_file(fR)
        dS = parse_file(fS)
        return np.array(dR), np.array(dS)

    def run(self):
        species=["AQIB01","CABD02","AQIA01","JYKQ01","ABDC01",
                 "ADFV01","AJFE01","AHZZ01","ABGA01","JZKE01",
                 "JABR01","AGCE01", "ABRT01"]
        loci= ['ighv','iglv', 'igkv','trav','trbv','trdv','trgv']

        dR, dS = self.get_data(fR, fS)


        #f, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        sns.set_palette("hls")
        sns.set_style("whitegrid")
        sns.distplot(dR, color="m") #, ax=axes[0, 0])
        sns.distplot(dS, color="g") #, ax=axes[0, 1])



        plt.tight_layout()
        plt.show()
# -----------------------------------------------
if __name__ == '__main__':

    #fR = "./trees_orig/res_probR.dat"
    #fS = "./trees_orig/res_probS.dat"


    fR = "./trees/res_probR.dat"
    fS = "./trees/res_probS.dat"

    L = LociRSdata()
    L.run()


