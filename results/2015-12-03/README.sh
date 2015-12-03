#!/bin/bash
#
#					2015-12-03
#					----------
#
# Here, I will run pyrad on the data. First I need to merge reads of the same sample from
# different sequencing runs.

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
              --min-overlap 12 \
              --threads 1 \
              --memory 800M &

         pear -f $DATADIR/../Coregonus2/${SAMPLE[$j]}_S$j'_L001_R1_001.fastq.gz' \
              -r $DATADIR/../Coregonus2/${SAMPLE[$j]}_S$j'_L001_R2_001.fastq.gz' \
              --output merged/${SAMPLE[$j]}_2 \
              --threads 1 \
              --memory 800M &
      done
      wait
   done
fi
