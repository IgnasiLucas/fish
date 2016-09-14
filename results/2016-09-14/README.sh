#!/bin/bash
#
#				2016-09-14
#				----------
#
# Here I just want to check how reliable the first hit of a merged pair of reads
# to the artificial contigs is. That is, how much better a merged pair of reads
# map to the first than to the second locus. I will use the index from 2016-06-29,
# and the reads 2016-05-03 (trimmed).


DATADIR=`pwd | sed 's/results/data/'`
FASTQDIR=`pwd | sed 's/2016-09-14/2016-05-03/'`
INDEXDIR=`pwd | sed 's/2016-09-14/2016-06-29/'`

if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi

SAMPLE=(ZeroNotUsed St0001 St0003        St0015 St0016 St0019 St0037 St0039
                    St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
                    Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

for i in `seq 1 23`; do
   if [ ! -e $DATADIR/${SAMPLE[$i]}.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed.fastq $DATADIR/${SAMPLE[$i]}.fastq
   fi
done

if [ ! -e fish.1.bt2 ]; then
   if [ -e $INDEXDIR/fish.1.bt2 ]; then
      ln -s $INDEXDIR/fish* ./
   else
      echo "Bowtie2 indices not found."
      exit
   fi
fi

function map_single_ends {
   bowtie2 --very-sensitive \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1 \
           --all \
           -x fish -U $DATADIR/$1.fastq \
           -S $1.se.sam &> $1.se.bowtie2.log
}

#for i in `seq 1 23`; do
for i in 13; do
   if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
      map_single_ends ${SAMPLE[$i]}
   fi
done

