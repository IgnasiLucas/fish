#!/bin/bash
#
#				2016-04-11
#				----------
# Given that a substantial portion of non-merged reads have an excess of Ns,
# it is possible that many of them would cluster with merged-reads loci. This
# is worth checking, since it would increase the coverage in those loci.
#
# If non-merged reads were coming, largely, from the same loci as merged
# reads, I would expect a negative correlation between the numbers of merged
# and non-merged reads, across samples. Or a relationship between the average
# quality and the proportion of non-merged reads. However, I don't see that.
# (Not shown).
#
# One option is to change the names of the reads in non-merged fastq files,
# to identify their origin. Then, I could use vsearch directly to check how
# many clusters are formed with merged and non-merged reads. As a simple test,
# I do not need to run all samples, but just one. They all have behaved quite
# similarly (except St0006), till now. I choose BlCl0091.
#

if [ ! -e BlCl0091.fasta ]; then
   cp ../2016-02-23/se/edits/BlCl0091.derep BlCl0091.fasta
   gawk '(/^>/){sub(/>/,">R1_"); print}(/^[^>]/){print}' ../2016-02-23/pe/edits/BlCl0091.firsts >> BlCl0091.fasta
fi

if [ ! -e BlCl0091.u ]; then
   vsearch --cluster_size BlCl0091.fasta \
           --threads 48 \
           --leftjust \
           --query_cov 50 \
           --id .86 \
           --userout BlCl0091.u \
           --userfields query+target+id+gaps+qstrand+qcov \
           --maxaccepts 1 \
           --maxrejects 0 \
           --minsl 0.2 \
           --usersort \
           --fulldp \
           --notmatched BlCl0091._temp
fi
