#!/bin/bash
#
#				2016-06-29
#				----------
#
# The consensus sequences of the read pairs that could be merged are ready in
# 2016-05-10/clust.86/pooled.consens.gz. It's a fasta file. Now, the idea is to
# use it as a reference genome, and map to it all the reads. Ideally, every
# consensus sequence should be covered by at least one merged read. The success
# of the approach can be measured by the proportion of non-merged reads that also
# map (unambiguously) to the consensus sequences, because these would be reads
# that failed to merge due to errors and could not be clustered with their
# homologous reads before.
#
# The fastq files, merged and trimmed, are in 2016-05-03.

DATADIR=`pwd | sed 's/results/data/'`
FASTQDIR=`pwd | sed 's/2016-06-29/2016-05-03/'`
if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi

if [ ! -e $DATADIR/consensus.fa ]; then
   cp ../2016-05-10/clust.86/pooled.consens.gz $DATADIR/
   gunzip $DATADIR/pooled.consens.gz
   mv     $DATADIR/pooled.consens $DATADIR/consensus.fa
fi

SAMPLE=(ZeroNotUsed St0001 St0003        St0015 St0016 St0019 St0037 St0039
                    St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
                    Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

for i in `seq 1 23`; do
   if [ ! -e $DATADIR/${SAMPLE[$i]}.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed.fastq $DATADIR/${SAMPLE[$i]}.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R1.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R1.fastq $DATADIR/${SAMPLE[$i]}_R1.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R2.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R2.fastq $DATADIR/${SAMPLE[$i]}_R2.fastq
   fi
done

if [ ! -e fish.1.bt2 ]; then
   bowtie2-build -seed 3817 --threads 4 $DATADIR/consensus.fa fish
fi

for i in `seq 1 23`; do
   if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
      bowtie2 --very-sensitive \
              --rg-id ${SAMPLE[$i]} \
              --rg "PL:ILLUMINA" \
              --rg "DT:2016" \
              --rg "SM:"${SAMPLE[$i]} \
              --un ${SAMPLE[$i]}.se.unmapped.fastq \
              -x fish -U ${SAMPLE[$i]}.fastq \
              -S ${SAMPLE[$i]}.se.sam &> ${SAMPLE[$i]}.se.bowtie2.log &

   fi
done
wait

for i in `seq 1 23`; do
   if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
      bowtie2 --sensitive \
              --maxins 600 \
              --rg-id ${SAMPLE[$i]} \
              --rg "PL:ILLUMINA" \
              --rg "DT:2016" \
              --rg "SM:"${SAMPLE[$i]} \
              --un-conc ${SAMPLE[$i]}.pe.unmapped_R%.fastq \
              -x fish -1 ${SAMPLE[$i]}_R1.fastq -2 ${SAMPLE[$i]}_R2.fastq \
              -S ${SAMPLE[$i]}.pe.sam &> ${SAMPLE[$i]}.pe.bowtie2.log &
   fi
done
