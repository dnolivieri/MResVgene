#!/bin/bash

#makeblastdb -in vconsensus.fasta -title vsCdb -dbtype prot -out ./vsCdb
makeblastdb -in vcons.fasta -title vsCdb -dbtype prot -out ./vsCdb
