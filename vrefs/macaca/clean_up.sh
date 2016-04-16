#!/bin/bash
# cat $i | sed '/^\s*$/d'

lbar="Pan_troglodytes Macaca_mulatta Gorilla_gorilla"

for i in $lbar; do
    cat $(ls $i*_sequence.fa) | sed '/^\s*$/d' > ${i}_all.fa
done
