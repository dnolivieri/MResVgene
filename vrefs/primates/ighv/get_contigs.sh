#!/bin/bash

animals="Aotus_nancymaae_JYKP01_ighv.fasta Callithrix_jacchus_JRUL01_ighv.fasta Chlorocebus_AQIB0_ighv_ighv.fasta Gorilla_gorilla_CABD02_ighv.fasta Macaca_fascicularis_AQIA01_ighv.fasta Macaca_mulatta_AANU01_ighv.fasta Macaca_nemestrina_JZLF01_ighv.fasta Mandrillus_leucophaeus_JYKQ01_ighv.fasta Microcebus_murinus_ABDC02_ighv.fasta Nomascus_leucogenys_ADFV01_ighv.fasta Pan_paniscus_AJFE01_ighv.fasta Pan_troglodytes_AACZ03_ighv.fasta Papio_anubis_AHZZ01_ighv.fasta Pongo_abelii_ABGA01_ighv.fasta Propithecus_coquereli_JZKE01_ighv.fasta Rhinopithecus_roxellana_JABR01_ighv.fasta Saimiri_AGCE01_ighv.fasta Tarsius_syrichta_ABRT02_ighv.fasta"


for i in $animals; do 
   echo $i
   grep \> $i  | awk -F"|" '{print $3}' | sort|uniq > ${i/.fasta/_cntg.txt}
done