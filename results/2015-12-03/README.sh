#!/bin/bash
#
#					2015-12-03
#					----------
#
# Here, I will merge the paired ends that can be merged, and pool reads from the same
# sample and alternative runs. I will get some statistics, and prepare the fastq files
# for pyrad.

DATADIR=`pwd | sed 's/results/data/'`
if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi

SAMPLE=(Undetermined St0001 St0003 St0006 St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

if [ ! -d fastqc ]; then
   mkdir fastqc
   mkdir fastqc/run1
   mkdir fastqc/run2
   fastqc -o fastqc/run1 -f fastq -t 6 $DATADIR/../Coregonus1/*_L001_R*_001.fastq.gz
   fastqc -o fastqc/run2 -f fastq -t 6 $DATADIR/../Coregonus2/*_L001_R*_001.fastq.gz
fi

if [ ! -d merged ]; then
   mkdir merged
   for i in 1 4 7 10 13 16 19 22; do
      for j in `seq $i $(( $i + 2 ))`; do
         pear -f $DATADIR/../Coregonus1/${SAMPLE[$j]}_S$j'_L001_R1_001.fastq.gz' \
              -r $DATADIR/../Coregonus1/${SAMPLE[$j]}_S$j'_L001_R2_001.fastq.gz' \
              --output merged/${SAMPLE[$j]}_1 \
              --min-overlap 10 \
              --quality-threshold 15 \
              --threads 1 \
              --memory 800M &

         pear -f $DATADIR/../Coregonus2/${SAMPLE[$j]}_S$j'_L001_R1_001.fastq.gz' \
              -r $DATADIR/../Coregonus2/${SAMPLE[$j]}_S$j'_L001_R2_001.fastq.gz' \
              --output merged/${SAMPLE[$j]}_2 \
              --min-overlap 10 \
              --quality-threshold 15 \
              --threads 1 \
              --memory 800M &
      done
      wait
   done
fi

# After merging, I want to visualize some statistics: the numbers of read pairs that were
# merged, not merged, and discarded; and also the length distributions of the merged reads.

if [ ! -e numbers.txt ] || find merged/ -type f -newer numbers.txt | grep -q fastq; then
   echo -e "#Sample\tRun\tMerged\tNot_merged\tDiscarded\tTotal" > numbers.txt
   echo -e "#Sample\tRun\tMerged\tNot_merged\tDiscarded\tTotal" > proportions.txt
   for i in `seq 1 24`; do
      M1=`wc -l merged/${SAMPLE[$i]}_1.assembled.fastq | gawk '{print $1/4}'`
      N1=`wc -l merged/${SAMPLE[$i]}_1.unassembled.forward.fastq | gawk '{print $1/4}'`
      D1=`wc -l merged/${SAMPLE[$i]}_1.discarded.fastq | gawk '{print $1/4}'`
      printf "%s\t%d\t%d\t%d\t%d\t%d\n" ${SAMPLE[$i]} 1 $M1 $N1 $D1 $(( $M1 + $N1 + $D1 )) >> numbers.txt

      M2=`wc -l merged/${SAMPLE[$i]}_2.assembled.fastq | gawk '{print $1/4}'`
      N2=`wc -l merged/${SAMPLE[$i]}_2.unassembled.forward.fastq | gawk '{print $1/4}'`
      D2=`wc -l merged/${SAMPLE[$i]}_2.discarded.fastq | gawk '{print $1/4}'`
      printf "%s\t%d\t%d\t%d\t%d\t%d\n" ${SAMPLE[$i]} 2 $M2 $N2 $D2 $(( $M2 + $N2 + $D2 )) >> numbers.txt
   done
fi

# I determined manually that the minimum and maximum lengths of merged reads are 50 and 592.
#for i in `seq 1 24`; do
#   if [ ! -e merged/${SAMPLE[$i]}_lengths.txt ] || find merged/ -type f -newer merged/${SAMPLE[$i]}_lengths.txt | grep -q fastq; then
#      gawk '(NR % 4 == 2){F[length($1)]++}END{for (i=50;i<=592;i++) print i "\t" F[i] + 0}' merged/${SAMPLE[$i]}_*.assembled.fastq > merged/${SAMPLE[$i]}_lengths.txt
#   fi
#done

if [ ! -e lengths.txt ] || find merged/ -type f -newer lengths.txt | grep -q .assembled.fastq; then
   HEADER="Length"
   seq 50 592 > zmatrix
   for i in `seq 1 24`; do
      gawk '(NR % 4 == 2){F[length($1)]++}END{for (i=50;i<=592;i++) print F[i] + 0}' merged/${SAMPLE[$i]}_*.assembled.fastq > zcolumn
      HEADER=$HEADER"\t"${SAMPLE[$i]}
      paste zmatrix zcolumn > zmatrixbig
      mv zmatrixbig zmatrix
   done
   echo -e $HEADER > zheader
   cat zheader zmatrix > lengths.txt
   rm zheader zmatrix zcolumn
fi

if [ ! -e merging.png ] || find . -type f -newer merging.png | grep -q numbers.txt; then
   R --save < plot01.R
fi

if [ ! -e sizes.png ] || find merged/ -type f -newer sizes.png | grep -q lengths.txt; then
   R --save < plot02.R
fi



# Conclusions
# -----------
#
# Merging was not very efficient. After changing the p-value of PEAR to 0.05 (more permissive
# than the default 0.01), the increase in proportion of merged read pairs was modest. I must
# work with both merged and unmerged reads.
#
# The size distribution is not identical among samples, but many samples have it similar,
# which is promising.
#
