#!/bin/bash


# ---- TRAV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep trav | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do    grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9]; done > Vextr_trav.dat
cat Vextr_trav.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_trav2.dat 

# --- TRBV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep trbv | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do  grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9];done > Vextr_trbv.dat
cat Vextr_trbv.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_trbv2.dat 


# --- TRDV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep trdv | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do  grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9];done > Vextr_trdv.dat
cat Vextr_trdv.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_trdv2.dat 

# --- IGHV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep ighv | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9];done > Vextr_ighv.dat
cat Vextr_ighv.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_ighv2.dat 


# --- IGKV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep igkv | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do    grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9]; done > Vextr_igkv.dat
cat Vextr_igkv.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_igkv2.dat 


# --- IGLV
lbar=$(grep \> Macaca_mulatta_MMUL01_outV_mrfq.fasta | grep iglv | awk -F"-" '{split($1,a,">"); print a[2]}')
for i in $lbar; do grep \> Macaca_mulatta_MMUL01_outVx.fasta | grep $i[^0-9]; done > Vextr_iglv.dat
cat Vextr_iglv.dat | awk -F"|" '{split($1,a,"["); print a[1],":"$2}' | awk -F">" '{print $2}' > Vextr_iglv2.dat 





