#!/bin/bash


grep \> Macaca_mulatta_MMUL1_ighv_outRF.fasta | awk -F"|" 'BEGIN{L=5000}
      {  split($1,a,"-"); gsub(/>/,"",a[1]);  
        gsub(/\(/, "",$2); gsub(/\)/, "",$2); 
        split($2,b,","); 
        if ($5==-1) {
          b[1]=L-b[1]; 
          b[2]=L-b[2]
        }; 
       print a[1], b[1], b[2], $5
     }' > vmres_ighv.dat



grep \> Macaca_mulatta_MMUL1_ighv_outV.fasta | awk -F"|"  '
  {
     a=$1
     gsub(/>/,"",a)
     print a, $2
  }'  > vext_ighv.dat
