#!/bin/bash
#
#				2015-12-16b
#				-----------
#
# I step back, starting from original fastq files, I use cutadapt to trim
# low quality tails, remove adapters, and filter out reads with unexpected
# barcode. Two surprises along the way: cutadapt cannot de-multiplex paired-
# end data, yet. And PEAR does not work with paired reads of different lengths.


DATADIR=`pwd | sed 's/results/data/'`
SAMPLE=(ZeroNotUsed St0001 St0003 St0006 St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)
BARCODE=(ZeroNotUsed AGCTA TCGAG TAACCTG ACCGAGT TGGGTGCC GTCTTGCG
         CAATCC GTTCAA AGCTA TCGAG TAACCTG ACCGAGT
         TGGGTGCC GTCTTGCG CAATCC GTTCAA AGCTA TCGAG
         TAACCTG ACCGAGT TGGGTGCC GTCTTGCG CAATCC GTTCAA)

if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi
if [ ! -d joined ]; then mkdir joined; fi
if [ ! -d merged ]; then mkdir merged; fi
if [ ! -d demultiplexed ]; then mkdir demultiplexed; fi
if [ ! -d trimmed ]; then mkdir trimmed; fi


# Joining of run 1 and run 2.
for j in 1 4 7 10 13 16 19 22; do
   for i in `seq $j $(( $j + 2 ))`; do
      if [ ! -e joined/${SAMPLE[$i]}_R1.fastq.gz ]; then
         cat_seqs -t fastq -z -o joined/${SAMPLE[$i]}_R1.fastq.gz \
                  $DATADIR/../Coregonus1/${SAMPLE[$i]}_S$i"_L001_R1_001.fastq.gz" \
                  $DATADIR/../Coregonus2/${SAMPLE[$i]}_S$i"_L001_R1_001.fastq.gz" &
      fi
      if [ ! -e joined/${SAMPLE[$i]}_R2.fastq.gz ]; then
         cat_seqs -t fastq -z -o joined/${SAMPLE[$i]}_R2.fastq.gz \
                  $DATADIR/../Coregonus1/${SAMPLE[$i]}_S$i"_L001_R2_001.fastq.gz" \
                  $DATADIR/../Coregonus2/${SAMPLE[$i]}_S$i"_L001_R2_001.fastq.gz" &
      fi
   done
   wait
done

# Merging with PEAR. From now on, I will have *.assembled.fastq and
# *.unassembled.forward.fastq and *.unassembled.reverse.fastq files
# for each sample. The '-q 5' option ends up cutting many forward reads
# on the first base.
for j in 1 7 13 19; do
   for i in `seq $j $(( $j + 5 ))`; do
      if [ ! -e merged/${SAMPLE[$i]}.assembled.fastq ]; then
         pear -f joined/${SAMPLE[$i]}_R1.fastq.gz \
              -r joined/${SAMPLE[$i]}_R2.fastq.gz \
              --output merged/${SAMPLE[$i]} \
              --min-overlap 12 \
              --threads 1 \
              --memory 800M > merged/${SAMPLE[$i]}_pear.log &
      fi
   done
   wait
done


# Unfortunately, cutadapt does not support de-multiplexing
# of paired-end files. I will demultiplex with SABRE first, and trim
# afterwards.
for j in 1 4 7 10 13 16 19 22; do
   for i in `seq $j $(( $j + 2 ))`; do
      if [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ]; then
         echo "${BARCODE[$i]} demultiplexed/${SAMPLE[$i]}_R1.fastq demultiplexed/${SAMPLE[$i]}_R2.fastq" > demultiplexed/${SAMPLE[$i]}.pe.barcode
         echo "${BARCODE[$i]} demultiplexed/${SAMPLE[$i]}.fastq" > demultiplexed/${SAMPLE[$i]}.barcode
         sabre pe -m 1 \
                  -f merged/${SAMPLE[$i]}.unassembled.forward.fastq \
                  -r merged/${SAMPLE[$i]}.unassembled.reverse.fastq \
                  -b demultiplexed/${SAMPLE[$i]}.pe.barcode \
                  -u demultiplexed/${SAMPLE[$i]}_unknown_R1.fastq \
                  -w demultiplexed/${SAMPLE[$i]}_unknown_R2.fastq > demultiplexed/${SAMPLE[$i]}_sabre.pe.log &

         sabre se -m 1 \
                  -f merged/${SAMPLE[$i]}.assembled.fastq \
                  -b demultiplexed/${SAMPLE[$i]}.barcode \
                  -u demultiplexed/${SAMPLE[$i]}_unknown.fastq > demultiplexed/${SAMPLE[$i]}_sabre.se.log &
      fi
   done
   wait
done

# Now, we trim assembled and unassembled reads separately as well. Hopefully,
# assembled reads should not need much trimming.
for j in 1 4 7 10 13 16 19 22; do
   for i in `seq $j $(( $j + 2 ))`; do
      if [ ! -e trimmed/${SAMPLE[$i]}_R1.fastq ]; then
         BARCODEREV=`echo ${BARCODE[$i]} | gawk -f revcomp.awk`
         cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                  -A "TA"$BARCODEREV"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
                  -o trimmed/${SAMPLE[$i]}_R1.fastq \
                  -p trimmed/${SAMPLE[$i]}_R2.fastq \
                  -q 17 \
                  demultiplexed/${SAMPLE[$i]}_R1.fastq \
                  demultiplexed/${SAMPLE[$i]}_R2.fastq > trimmed/${SAMPLE[$i]}_cutadapt.pe.log &

         cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                  -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT${BARCODE[$i]} \
                  -o trimmed/${SAMPLE[$i]}.fastq \
                  -q 17 \
                  demultiplexed/${SAMPLE[$i]}.fastq > trimmed/${SAMPLE[$i]}_cutadapt.se.log &
      fi
   done
   wait
done

## Clean up.
#for i in `seq 1 24`; do
#   if [ -e merged/${SAMPLE[$i]}.assembled.fastq ]; then
#      rm demultiplexed/${SAMPLE[$i]}*.fastq
#      rm trimmed/${SAMPLE[$i]}*.fastq.gz
#   fi
#done

cd ../2015-12-17
./README.sh 1> readme.log 2> readme.err
