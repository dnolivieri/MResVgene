#!/bin/bash
#
#  dnolivieri (18Feb2013)
#    - automating the PhyML 

#   This needed a fix for the writing of the phy file...
#  head -4 elephant_loci.phy | sed -e '2,4 !d' -e 's/\([^ ]\{10\}\)\([^ ].*\)/\1  \2/'

#   cat elephant_loci.phy | awk '{s=$0; if (s~ /^Vs[0-9]/){a=substr(s,1,10); b=substr(s,11); c=sprintf("%s  %s", a, b); print c}else{print s }}' 


for i in $(ls mammals/$1/$1*.fa ); do 
   echo $i 

   guideOut=${i/fa/dnd} 
   outfile=${i/fa/phy}

   ## For the study, this has been done... For other animals, uncomment.
   #cat $i | sed -e 's/\[/-/' -e 's/\]:/-/' > tmp.fa
   cat $i | awk '$0~">"{split($0,a,"|"); s=sprintf("%s-%s-%s",a[1],substr(a[2],1,6),substr(a[4],2,5)); print s }$0!~">"{print $0}' > tmp.fa



   mv tmp.fa $i

   filename=$i
   #echo $i
   #clustalo  --infile=$filename --guidetree-out=$guideOut --outfile=$outfile  --outfmt=phy -v  --force

   ### Note, my version of clustalo fixes the size of the name
   #/home/david/Research/BioInf/src/clustal-omega-1.1.0/src/clustalo  --infile=$filename --guidetree-out=$guideOut   --outfmt=phy --force | awk '{s=$0; if (s~ /^Vs[0-9]/){a=substr(s,1,10); b=substr(s,11); c=sprintf("%s  %s", a, b); print c}else{print s }}'  > $outfile

   /home/david/Research/BioInf/src/clustal-omega-1.1.0/src/clustalo  --infile=$filename --guidetree-out=$guideOut --outfile=$outfile  --outfmt=phy --force 

   phyml_infile=$outfile
   /home/david/Research/BioInf/seaview/PhyML_3.0_linux32 -d aa -m LG -b -4 -v 0.0 -c 4 -a e -f m -i $phyml_infile
   #sleep 1  ## I am wondering if this will help.
   #nohup blastp -out ${i/.fa/.xml}  -outfmt 5 -query $i -db /usr/local/ncbiNRdb/nr &
done

