#!/bin/bash

if [ ! -e z1 ]; then
   velveth assembly/ 21 -fastq -shortPaired Bl0091_trimmed_R1.fastq Bl0091_trimmed_R2.fastq -separate -strand_specific
   echo "velveth run" > z1
fi

if [ ! -e z2 ]; then
   velvetg assembly/ -cov_cutoff 1 -min_contig_lgth 75 -ins_length 500 -exp_cov 2 -scaffolding no -read_trkg yes -unused_reads yes
   echo "velvetg run" > z2
fi
