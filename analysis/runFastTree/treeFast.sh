#!/bin/bash
#
#  dnolivieri (17Abr2013)
#    - automating the FastTree 
#      used as a script in the web.
#To infer a tree for a protein alignment with the JTT+CAT model, use
#FastTree < alignment_file > tree_file 
#or
#FastTree alignment.file > tree_file 
#Use the -wag option to use the WAG+CAT model instead of JTT+CAT.
#To infer a tree for a nucleotide alignment with the GTR+CAT model, use
#FastTree -gtr -nt < alignment.file > tree_file 
#or
#FastTree -gtr -nt alignment_file > tree_file 
#If you do not specify -gtr, then FastTree will use the Jukes-Cantor + CAT model instead

#Does FastTree require aligned sequences?
#Yes. If you do not have a multiple sequence alignment, we recommend using MUSCLE to create one 
#and gblocks to "trim" it (remove low-confidence parts of the alignment, such as positions that 
#contain many gaps). Alternatively, for large protein families, we recommend hmmalign from the 
#HMMer package. If using hmmalign, you should remove columns that do not match the model (usually 
#output in lower case). For large RNA families we recommend Infernal.


  

  i=$1
  guideOut=${i/.fasta/.dnd} 
  outfile=${i/.fasta/.fasta}
  ft_outfile=${i/.fasta/.nwk}

<<EOF
  ## For the study, this has been done... For other animals, uncomment.
  #cat $i | sed -e 's/\[/-/' -e 's/\]:/-/' > /var/www/VgeneDB/site/web/phylo/fasttree/tmp.fa

  cat $i |  awk '{if ($0~/>/){ split($0,a,"|"); split(a[2],v,"_"); s=sprintf("%s-%s%s-%s",a[1], substr(v[1],1,4), toupper(substr(v[2],1,1)), substr(a[4],2,5)); print s}else{print $0}}'> tmp.fa

  mv tmp.fa $i
EOF

   filename=$i
   /media/disk2TB/BioInf/src/clustal-omega-1.1.0/src/clustalo  --infile=$filename --guidetree-out=$guideOut --outfile=$outfile  --outfmt=fasta --force 
   ./fasttree/FastTree -wag $outfile > $ft_outfile

