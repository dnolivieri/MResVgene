#!/bin/bash
##  run as ./runTrees.sh <Fast|ML> 


#Species="../primates_res/all_init.fasta"


#Species="Chlorocebus_AQIB0.fasta Gorilla_gorilla_CABD02.fasta Macaca_fascicularis_AQIA01.fasta Mandrillus_leucophaeus_JYKQ01.fasta Microcebus_murinus_ABDC01.fasta Nomascus_leucogenys_ADFV01.fasta Pan_paniscus_AJFE01.fasta Papio_anubis_AHZZ01.fasta Pongo_abelii_ABGA01.fasta Propithecus_coquereli_JZKE01.fasta Rhinopithecus_roxellana_JABR01.fasta Saimiri_AGCE01.fasta"
#Species="Chlorocebus_AQIB01.fasta"
#Species="Gorilla_gorilla_CABD02.fasta"
#Species="Macaca_fascicularis_AQIA01.fasta"


#Species="Tarsius_syrichta_ABRT01.fasta"
#srcdir1="/media/disk2TB/BioInf/VsRF/analysis/treeMacacaM"
#Species="Macaca_mulatta_MMUL01.fasta"
#Species="all_R.fasta all_S.fasta"
#srcdir1="/media/disk2TB/BioInf/VsRF/analysis/MacacaTrees"
#Species="macaca_trav.fasta  macaca_trbv.fasta"


#Species="Chlorocebus_AQIB01.fasta  Gorilla_gorilla_CABD02.fasta  Macaca_fascicularis_AQIA01.fasta Mandrillus_leucophaeus_JYKQ01.fasta Microcebus_murinus_ABDC01.fasta Nomascus_leucogenys_ADFV01.fasta Pan_paniscus_AJFE01.fasta Papio_anubis_AHZZ01.fasta Pongo_abelii_ABGA01.fasta Propithecus_coquereli_JZKE01.fasta Rhinopithecus_roxellana_JABR01.fasta Saimiri_AGCE01.fasta Tarsius_syrichta_ABRT01.fasta"
srcdir1="/media/disk2TB/BioInf/VsRF/analysis/cmpLociTrees"
Species="all_R.fasta all_S.fasta"
Species="AQIB01_Ighv.fasta AQIB01_igkv.fasta AQIB01_iglv.fasta AQIB01_trav.fasta AQIB01_trbv.fasta AQIB01_trgv.fasta"

for i in $Species; do 
   echo $i 
   sleep 1
   echo "fasttree";
   #nohup ./treeFast.sh $i &
   ./treeFast.sh $srcdir1/$i
done

 



