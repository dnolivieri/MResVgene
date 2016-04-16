#!/bin/bash


srcdir1=/media/disk2TB/BioInf/VsRF/bstrap
srcdir2=/media/disk2TB/BioInf/VsRF/vrefs/primates
f1="3_outRF.fasta"
f2="outV_mr.fasta"

mlist="Chlorocebus_AQIB01 Gorilla_gorilla_CABD02 Macaca_fascicularis_AQIA01 Mandrillus_leucophaeus_JYKQ01 Microcebus_murinus_ABDC01 Nomascus_leucogenys_ADFV01 Pan_paniscus_AJFE01 Papio_anubis_AHZZ01 Pongo_abelii_ABGA01 Propithecus_coquereli_JZKE01 Rhinopithecus_roxellana_JABR01 Saimiri_AGCE01"

#mlist="Chlorocebus_AQIB01" 
#mlist="Gorilla_gorilla_CABD02"
#mlist="Macaca_fascicularis_AQIA01"
#mlist="Microcebus_murinus_ABDC01"

for m in $mlist; do 
    cat $srcdir1/${m}_$f1  $srcdir2/${m}_$f2 > ${m}.fasta
done; 




















