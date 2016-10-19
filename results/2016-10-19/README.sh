#!/bin/bash
#
#				2016-10-19
#				----------
#
# The consensus sequences of the read pairs that could be merged are ready in
# 2016-10-03/clust.86/pooled.consens.gz. It's a fasta file. Now, the idea is to
# use it as a reference genome, and map to it all the reads. Ideally, every
# consensus sequence should be covered by at least one merged read. The success
# of the approach can be measured by the proportion of non-merged reads that also
# map (unambiguously) to the consensus sequences, because these would be reads
# that failed to merge due to errors and could not be clustered with their
# homologous reads before.
#
# The fastq files are in 2016-10-03 (trimmed paired ends) or in 2016-10-03/merged
# (merged reads, non-trimmed).

PREVIOUS=`pwd | sed 's/2016-10-19/2016-06-29/'`
 DATADIR=`pwd | sed 's/results/data/'`
FASTQDIR=`pwd | sed 's/2016-10-19/2016-10-03/'`
if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi

if [ ! -e $DATADIR/consensus.fa ]; then
   cp ../2016-10-03/clust.86/pooled.consens.gz $DATADIR/
   gunzip $DATADIR/pooled.consens.gz
   mv     $DATADIR/pooled.consens $DATADIR/consensus.fa
fi

SAMPLE=(ZeroNotUsed StCa0001 StCa0003          StCa0015 StCa0016 StCa0019 StCf0037 StCf0039
                    StCf0043 StCf0044 StCf0049 StCf0050 BlCa0065 BlCa0076 BlCa0080 BlCa0083
                    BlCa0104 BlCa0108 BlCl0091 BlCl0093 BlCl0094 BlCl0095 BlCl0098 BlCl0116)

for i in `seq 1 23`; do
   if [ ! -e $DATADIR/${SAMPLE[$i]}.fastq ]; then
      ln -s $FASTQDIR/merged/${SAMPLE[$i]}.assembled.fastq $DATADIR/${SAMPLE[$i]}.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R1.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R1.fastq $DATADIR/${SAMPLE[$i]}_R1.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R2.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R2.fastq $DATADIR/${SAMPLE[$i]}_R2.fastq
   fi
done

if [ ! -e fish.1.bt2 ]; then
   bowtie2-build $DATADIR/consensus.fa fish
fi

function map_single_ends {
   bowtie2 --sensitive \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1 \
           --un $1.se.unmapped.fastq \
           -x fish -U $DATADIR/$1.fastq \
           -S $1.se.sam &> $1.se.bowtie2.log
}

function map_paired_ends {
   bowtie2 --sensitive \
           --maxins 600 \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1} \
           --un-conc $1.pe.unmapped_R%.fastq \
           -x fish -1 $DATADIR/$1'_R1.fastq' -2 $DATADIR/$1'_R2.fastq' \
           -S $1.pe.sam &> $1.pe.bowtie2.log
}

PROC=`grep -P '^processor' /proc/cpuinfo | wc -l`
#              23
if [ $PROC -gt 46 ]; then
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]} &
      fi
   done
#   wait
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]} &
      fi
   done
   wait
elif [ $PROC -ge 6 ]; then
   for j in 1 7 13; do
      for i in `seq $j $(( j + 5 ))`; do
         if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
            map_single_ends ${SAMPLE[$i]} &
         fi
      done
      wait
      for i in `seq $j $(( j + 5 ))`; do
         if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
            map_paired_ends ${SAMPLE[$i]} &
         fi
      done
      wait
   done
   for i in 19 20 21 22 23; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]} &
      fi
   done
   wait
   for i in 19 20 21 22 23; do
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]} &
      fi
   done
   wait
else
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]}
      fi
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]}
      fi
   done
fi

#if [ ! -e summary_se.txt ]; then
   echo "## Tot.Reads: Total number of reads."                             > summary_se.txt
   echo "## Unmapped:  Aligned 0 times."                                  >> summary_se.txt
   echo "## Uniq.Map.: Aligned exactly 1 time."                           >> summary_se.txt
   echo "## Mult.Map.: Aligned >1 times."                                 >> summary_se.txt
   echo "## Overall: Overall mapping rate."                               >> summary_se.txt
   echo -e "# Sample\tTot.Reads\tUnmapped\tUniq.Map.\tMult.Map.\tOverall" >> summary_se.txt
   gawk -f summary_se.awk *.se.bowtie2.log | sort -nrk 2                  >> summary_se.txt
#fi
#if [ ! -e summary_pe.txt ]; then
   echo "## Tot.Pairs: Total number of read pairs."                        > summary_pe.txt
   echo "## Conc.Uniq.: Uniquely and concordantly mapped pairs."          >> summary_pe.txt
   echo "## Conc.Mult.: Multiply and concordantly mapped pairs."          >> summary_pe.txt
   echo "## Discord.:  Discordantly mapped pairs."                        >> summary_pe.txt
   echo "## Un.Pairs: Unmapped pairs."                                    >> summary_pe.txt
   echo "## Un.Mates: Reads or mates from unmapped pairs that do not map either as single reads." >> summary_pe.txt
   echo "## Uniq.Mates: Mates from unmapped pairs that map uniquely as single reads."             >> summary_pe.txt
   echo "## Mult.Mates: Mates from unmapped pairs that map ambiguously as single reads."          >> summary_pe.txt
   echo "## Overall: Overall mapping rate."                               >> summary_pe.txt
   echo -e "# Sample\tTot.Pairs\tConc.Uniq.\tConc.Mult.\tDiscord.\tUn.Pairs\tUn.Mates\tUniq.Mates\tMult.Mates\tOverall" >> summary_pe.txt
   gawk -f summary_pe.awk *.pe.bowtie2.log | sort -nrk 2                  >> summary_pe.txt
#fi
gnuplot < summary_se.gnuplot
gnuplot < summary_pe.gnuplot

# This is to compare the current results with those from 2016-06-29. Because
# samples behave similarly between the two clustering experiments, I compare
# only the total numbers of reads.

#if [ ! -e comparison.txt ]; then
   echo "                       "  > z1
   echo "Total merged reads     " >> z1
   echo "   Uniquely mapped     " >> z1
   echo "   Ambiguously mapped  " >> z1
   echo "   Unmapped            " >> z1
   echo                           >> z1
   echo "Total pairs            " >> z1
   echo "   Concordant, unique  " >> z1
   echo "   Concor., ambiguous  " >> z1
   echo "   Discordant          " >> z1
   echo "   Unmapped            " >> z1
   echo                           >> z1
   echo "Mates in unmapped pairs" >> z1
   echo "   Uniquely mapped     " >> z1
   echo "   Ambiguously mapped  " >> z1
   echo "   Unmapped            " >> z1

   basename $PREVIOUS > z2
   gawk -f totals.awk $PREVIOUS/summary*.txt >> z2

   basename `pwd` > z3
   gawk -f totals.awk ./summary*.txt >> z3

   paste z1 z2 z3 > comparison.txt
   rm z1 z2 z3
#fi

# -------------------------------------------------------------------------
#                       	     2016-06-29	           2016-10-19
# -------------------------------------------------------------------------
# Total merged reads     	12221477          	12221553
#    Uniquely mapped     	 5746594 ( 47.02 %)	 5641845 ( 46.16 %)
#    Ambiguously mapped  	 5076638 ( 41.54 %)	 3321613 ( 27.18 %)
#    Unmapped            	 1398245 ( 11.44 %)	 3258095 ( 26.66 %)
#
# Total pairs            	15671492          	15671492
#    Concordant, unique  	 2694352 ( 17.19 %)	 1542739 (  9.84 %)
#    Concor., ambiguous  	  804347 (  5.13 %)	  198427 (  1.27 %)
#    Discordant          	  218040 (  1.39 %)	  408478 (  2.61 %)
#    Unmapped            	11954753 ( 76.28 %)	13521848 ( 86.28 %)
#
# Mates in unmapped pairs	23909506          	27043696
#    Uniquely mapped     	 2865898 ( 11.99 %)	 3951440 ( 14.61 %)
#    Ambiguously mapped  	10004087 ( 41.84 %)	 8694552 ( 32.15 %)
#    Unmapped            	11039521 ( 46.17 %)	14397704 ( 53.24 %)
# -------------------------------------------------------------------------
#
#
# CONCLUSIONS
# -----------
#
# The current adjustment of parameters somewhat improved the clustering of
# merged, pooled reads, but did not get rid of the split loci, as can be seen
# from the still large portion of ambiguously mapped merged reads (27.18 %).
# Therefore, it is urgent to move on along the following terms:
#
#    1. Give up the clusters of pooled merged reads as a potential reference
#       to identify variants.
#
#    2. Discard any paired-ends that map in any way to the merged reads in
#       the first analysis, and analyse the remaining 76% of unmapped paired
#       ends separately from the merged reads, using ipyrad.
#
