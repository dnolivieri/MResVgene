#!/bin/bash


#srcdir1=/media/disk2TB/BioInf/VsRF/bstrap
#srcdir1=/media/disk2TB/BioInf/VsFW/analysis_bstrap
srcdir1=/media/disk2TB/BioInf/VsFW/analysis_bstrapNew/iter4
srcdir2=/media/disk2TB/BioInf/VsRF/vrefs/primates
f1="4_outRF_mrfq.fasta"
f2="outV_mrfq.fasta"

mlist="Chlorocebus_AQIB01 Gorilla_gorilla_CABD02 Macaca_fascicularis_AQIA01 Mandrillus_leucophaeus_JYKQ01 Microcebus_murinus_ABDC01 Nomascus_leucogenys_ADFV01 Pan_paniscus_AJFE01 Papio_anubis_AHZZ01 Pongo_abelii_ABGA01 Propithecus_coquereli_JZKE01 Rhinopithecus_roxellana_JABR01 Saimiri_AGCE01 Tarsius_syrichta_ABRT01"

#mlist="Chlorocebus_AQIB01" 
#mlist="Gorilla_gorilla_CABD02"
#mlist="Macaca_fascicularis_AQIA01"
#mlist="Microcebus_murinus_ABDC01"

for m in $mlist; do 
    cat $srcdir1/${m}_$f1  $srcdir2/${m}_$f2 > ${m}.fasta
done; 




















