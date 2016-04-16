#!/bin/bash


lbar="Chinchilla_lanigera_AGCD01_outV.fasta Cricetulus_griseus_AFTD01_outV.fasta Dipodomys_ordii_ABRO01_outV.fasta Fukomys_damarensis_AYUG01_outV.fasta Heterocephalus_glaber_AHKG01_outV.fasta Ictidomys_tridecemlineatus_AGTP01_outV.fasta Jaculus_jaculus_AKZC01_outV.fasta Mesocricetus_auratus_APMT01_outV.fasta Microtus_ochrogaster_AHZW01_outV.fasta Nannospalax_galili_AXCS01_outV.fasta Octodon_degus_AJSA01_outV.fasta Peromyscus_maniculatus_AYHN01_outV.fasta"

for i in $lbar; do 
    echo $i ${i/.fasta/_cntg.txt}
    grep \> $i | awk -F"|" '{print $3}' |sort |uniq > ${i/.fasta/_cntg.txt}
done;