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
           --threads 64 \
           --leftjust \
           --id 0.86 \
           --userout BlCl0091.u \
           --userfields query+target+id+gaps+qstrand+qcov \
           --maxaccepts 1 \
           --maxrejects 0 \
           --minsl 0.2 \
           --notmatched BlCl0091._temp
fi

if [ ! -e summary.txt ]; then
   gawk 'BEGIN{
      TRANSLATE["R1Bl"] = "Cross-hit          "
      TRANSLATE["BlR1"] = "Cross-hit          "
      TRANSLATE["R1R1"] = "Non-merged self-hit"
      TRANSLATE["BlBl"] = "Merged self-hit    "
      print               "# Type_of_hit      \tFreq.\tIdent.\tGaps"
   }{
      A = substr($1,1,2)
      B = substr($2,1,2)
      F[TRANSLATE[A B]]++
      IDENTITY[TRANSLATE[A B]] += $3
      GAPS[TRANSLATE[A B]] += $4
   }END{
      for (f  in F) {
         print f "\t" F[f] "\t" IDENTITY[f]/F[f] "\t" GAPS[f]/F[f]
      }
   }' BlCl0091.u > summary.txt

   R1SINGLE=`grep ">R1" BlCl0091._temp | wc -l`
   BLSINGLE=`grep ">Bl" BlCl0091._temp | wc -l`
   echo -e "Non-merged singletons\t$R1SINGLE\t-\t-" >> summary.txt
   echo -e "Merged singletons\t$BLSINGLE\t-\t-" >> summary.txt
fi

# Conclusion
# ----------
#
# ---------------------------------------------
# Type_of_hit      	Freq.	Ident.	Gaps
# ---------------------------------------------
# Non-merged self-hit	261056	92.9803	2.43714
# Merged self-hit    	202574	94.223	4.4491
# Cross-hit          	219029	92.1468	5.8409
# Non-merged singletons	209932	-	-
# Merged singletons	281981	-	-
# ---------------------------------------------
#
# Among the 951073 non-merged (unique) reads from BlCl0091, 219029 (23%)
# have their most similar read among merged reads. This is an important
# number of reads that can increase the coverage of several loci. However,
# the number of non-merged-specific loci is even higher, and should not
# be discarded.

