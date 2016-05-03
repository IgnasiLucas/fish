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
      if [ ! -e merged/${SAMPLE[$i]}.assembled.fastq ] && [ ! -e merged/${SAMPLE[$i]}.assembled.fastq.gz ]; then
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
# Clean up. I remove discarded reads.
find merged -name '*.discarded.fastq' -exec rm '{}' \;

# 2016-04-13: I just learned that PEAR has reverse-complemented second reads.
# This has implications in downstream analysis (trimming of adapters, and pyrad).
# I think it is better to get the second reads back in their original orientation.
touch checkpoints
for i in `seq 1 24`; do
   if ! grep -q "Reverse ${SAMPLE[$i]}" checkpoints; then
      echo "Reverse ${SAMPLE[$i]}" >> checkpoints
      # Executing this in the server, I want each pair of commands to be run in
      # parallel. That is, mv should not be executed before vsearch finishes, within
      # each sample, but all samples can be run together. Thus, I use curly brackets.
      {
         vsearch --fastx_revcomp merged/${SAMPLE[$i]}.unassembled.reverse.fastq \
                 --fastqout merged/${SAMPLE[$i]}.unassembled.reverse.reversed.fastq
         mv merged/${SAMPLE[$i]}.unassembled.reverse.reversed.fastq merged/${SAMPLE[$i]}.unassembled.reverse.fastq
      } &
   fi
done
wait

## Unfortunately, cutadapt does not support de-multiplexing
## of paired-end files. I will demultiplex with SABRE first, and trim
## afterwards.
#for j in 1 4 7 10 13 16 19 22; do
#   for i in `seq $j $(( $j + 2 ))`; do
#      if [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
#         echo "${BARCODE[$i]} demultiplexed/${SAMPLE[$i]}_R1.fastq demultiplexed/${SAMPLE[$i]}_R2.fastq" > demultiplexed/${SAMPLE[$i]}.pe.barcode
#         echo "${BARCODE[$i]} demultiplexed/${SAMPLE[$i]}.fastq" > demultiplexed/${SAMPLE[$i]}.barcode
#         sabre pe -m 1 \
#                  -f merged/${SAMPLE[$i]}.unassembled.forward.fastq \
#                  -r merged/${SAMPLE[$i]}.unassembled.reverse.fastq \
#                  -b demultiplexed/${SAMPLE[$i]}.pe.barcode \
#                  -u demultiplexed/${SAMPLE[$i]}_unknown_R1.fastq \
#                  -w demultiplexed/${SAMPLE[$i]}_unknown_R2.fastq > demultiplexed/${SAMPLE[$i]}_sabre.pe.log &
#
#         sabre se -m 1 \
#                  -f merged/${SAMPLE[$i]}.assembled.fastq \
#                  -b demultiplexed/${SAMPLE[$i]}.barcode \
#                  -u demultiplexed/${SAMPLE[$i]}_unknown.fastq > demultiplexed/${SAMPLE[$i]}_sabre.se.log &
#      fi
#   done
#   wait
#done

# 2016-04-13: Before, I used Sabre. It works, but now that I have to repeat the
# demultiplexing, after having reversed-complemented the non-merged second reads
# back to their original orientation, I will use pyrad for individual sample
# demultiplexing, because it gives more information in the stats. Note that now
# I run this in the server, all samples at once.
for i in `seq 1 24`; do
   if [ ! -d demultiplexed/${SAMPLE[$i]} ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ]; then mkdir demultiplexed/${SAMPLE[$i]}; fi
   if [ ! -d demultiplexed/${SAMPLE[$i]}/pe_full ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ]; then mkdir demultiplexed/${SAMPLE[$i]}/pe_full; fi
   if [ ! -d demultiplexed/${SAMPLE[$i]}/pe_trunc ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ]; then mkdir demultiplexed/${SAMPLE[$i]}/pe_trunc; fi
   if [ ! -d demultiplexed/${SAMPLE[$i]}/se_full ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ]; then mkdir demultiplexed/${SAMPLE[$i]}/se_full; fi
   if [ ! -d demultiplexed/${SAMPLE[$i]}/se_trunc ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ]; then mkdir demultiplexed/${SAMPLE[$i]}/se_trunc; fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/barcode.1 ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ]; then
      echo -e ${SAMPLE[$i]}"\t"${BARCODE[$i]} > demultiplexed/${SAMPLE[$i]}/barcode.1
   fi
   # I notice that many reads miss the first base, and the truncated barcode is not
   # recognized. I will run pyrad's step 1 twice, once with the complete barcode and
   # 1 mismatch allowed, and once with the truncated barcode, without any mismatch allowed.
   if [ ! -e demultiplexed/${SAMPLE[$i]}/barcode.2 ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      gawk '{print $1 "\t" substr($2,2)}' demultiplexed/${SAMPLE[$i]}/barcode.1 > demultiplexed/${SAMPLE[$i]}/barcode.2
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R1_.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      ln -s `pwd`"/merged/${SAMPLE[$i]}.unassembled.forward.fastq" demultiplexed/${SAMPLE[$i]}_R1_.fastq
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R2_.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R2.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R2.fastq.gz ]; then
      ln -s `pwd`"/merged/${SAMPLE[$i]}.unassembled.reverse.fastq" demultiplexed/${SAMPLE[$i]}_R2_.fastq
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/params_1.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      pyrad -n
      sed -i '/## 1. /c\demultiplexed/'${SAMPLE[$i]}'/pe_full/   ## 1. Working directory                                 (all)' params.txt
      sed -i '/## 2. /c\demultiplexed/'${SAMPLE[$i]}'*.fastq     ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)' params.txt
      sed -i '/## 3. /c\demultiplexed/'${SAMPLE[$i]}'/barcode.1  ## 3. Loc. of barcode file (if not line 18)             (s1)' params.txt
      sed -i '/## 6. /c\TAA,CGG                                  ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\1                                        ## 7. N processors (parallel)                           (all)' params.txt
      sed -i '/## 11./c\pairddrad                                ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)' params.txt
      sed -i '/## 19./c\1                                        ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)' params.txt
      mv params.txt demultiplexed/${SAMPLE[$i]}/params_1.txt
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/params_2.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      pyrad -n
      sed -i '/## 1. /c\demultiplexed/'${SAMPLE[$i]}'/pe_trunc/  ## 1. Working directory                                 (all)' params.txt
      sed -i '/## 2. /c\demultiplexed/'${SAMPLE[$i]}'*.fastq     ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)' params.txt
      sed -i '/## 3. /c\demultiplexed/'${SAMPLE[$i]}'/barcode.2  ## 3. Loc. of barcode file (if not line 18)             (s1)' params.txt
      sed -i '/## 6. /c\TAA,CGG                                  ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\1                                        ## 7. N processors (parallel)                           (all)' params.txt
      sed -i '/## 11./c\pairddrad                                ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)' params.txt
      sed -i '/## 19./c\0                                        ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)' params.txt
      mv params.txt demultiplexed/${SAMPLE[$i]}/params_2.txt
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/params_3.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      pyrad -n
      sed -i '/## 1. /c\demultiplexed/'${SAMPLE[$i]}'/se_full/   ## 1. Working directory                                 (all)' params.txt
      sed -i '/## 2. /c\merged/'${SAMPLE[$i]}'.assembled.fastq   ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)' params.txt
      sed -i '/## 3. /c\demultiplexed/'${SAMPLE[$i]}'/barcode.1  ## 3. Loc. of barcode file (if not line 18)             (s1)' params.txt
      sed -i '/## 6. /c\TAA,CGG                                  ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\1                                        ## 7. N processors (parallel)                           (all)' params.txt
      sed -i '/## 11./c\ddrad                                    ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)' params.txt
      sed -i '/## 19./c\1                                        ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)' params.txt
      mv params.txt demultiplexed/${SAMPLE[$i]}/params_3.txt
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/params_4.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      pyrad -n
      sed -i '/## 1. /c\demultiplexed/'${SAMPLE[$i]}'/se_trunc/  ## 1. Working directory                                 (all)' params.txt
      sed -i '/## 2. /c\merged/'${SAMPLE[$i]}'.assembled.fastq   ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)' params.txt
      sed -i '/## 3. /c\demultiplexed/'${SAMPLE[$i]}'/barcode.2  ## 3. Loc. of barcode file (if not line 18)             (s1)' params.txt
      sed -i '/## 6. /c\TAA,CGG                                  ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\1                                        ## 7. N processors (parallel)                           (all)' params.txt
      sed -i '/## 11./c\ddrad                                    ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)' params.txt
      sed -i '/## 19./c\0                                        ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)' params.txt
      mv params.txt demultiplexed/${SAMPLE[$i]}/params_4.txt
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/pe_full/stats/s1.sorting.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      pyrad -p demultiplexed/${SAMPLE[$i]}/params_1.txt -s 1 &
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/pe_trunc/stats/s1.sorting.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      pyrad -p demultiplexed/${SAMPLE[$i]}/params_2.txt -s 1 &
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/se_full/stats/s1.sorting.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      pyrad -p demultiplexed/${SAMPLE[$i]}/params_3.txt -s 1 &
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}/se_trunc/stats/s1.sorting.txt ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      pyrad -p demultiplexed/${SAMPLE[$i]}/params_4.txt -s 1 &
   fi
done
wait
for i in `seq 1 24`; do
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ]; then
      if [ `head -n 2 demultiplexed/${SAMPLE[$i]}/pe_trunc/stats/s1.sorting.txt | tail -n 1 | cut -f 4` -gt 0 ]; then
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}_R1.fastq \
            demultiplexed/${SAMPLE[$i]}/pe_full/fastq/${SAMPLE[$i]}_R1.fq.gz \
            demultiplexed/${SAMPLE[$i]}/pe_trunc/fastq/${SAMPLE[$i]}_R1.fq.gz &
      else
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}_R1.fastq \
            demultiplexed/${SAMPLE[$i]}/pe_full/fastq/${SAMPLE[$i]}_R1.fq.gz &
      fi
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R2.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}_R2.fastq.gz ]; then
      if [ `head -n 2 demultiplexed/${SAMPLE[$i]}/pe_trunc/stats/s1.sorting.txt | tail -n 1 | cut -f 4` -gt 0 ]; then
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}_R2.fastq \
            demultiplexed/${SAMPLE[$i]}/pe_full/fastq/${SAMPLE[$i]}_R2.fq.gz \
            demultiplexed/${SAMPLE[$i]}/pe_trunc/fastq/${SAMPLE[$i]}_R2.fq.gz &
      else
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}_R2.fastq \
            demultiplexed/${SAMPLE[$i]}/pe_full/fastq/${SAMPLE[$i]}_R2.fq.gz &
      fi
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ]; then
      if [ `head -n 2 demultiplexed/${SAMPLE[$i]}/se_trunc/stats/s1.sorting.txt | tail -n 1 | cut -f 4` -gt 0 ]; then
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}.fastq \
            demultiplexed/${SAMPLE[$i]}/se_full/fastq/${SAMPLE[$i]}_R1.fq.gz \
            demultiplexed/${SAMPLE[$i]}/se_trunc/fastq/${SAMPLE[$i]}_R1.fq.gz &
      else
         cat_seqs -t fastq -o demultiplexed/${SAMPLE[$i]}.fastq \
            demultiplexed/${SAMPLE[$i]}/se_full/fastq/${SAMPLE[$i]}_R1.fq.gz &
      fi
   fi
done
wait
# Clean up.
for i in `seq 1 24`; do
   if [ -e demultiplexed/${SAMPLE[$i]}_R1.fastq ] && [ -e demultiplexed/${SAMPLE[$i]}_R2.fastq ] && [ -e demultiplexed/${SAMPLE[$i]}.fastq ] && [ -d demultiplexed/${SAMPLE[$i]} ]; then
      rm -r demultiplexed/${SAMPLE[$i]}
      find demultiplexed -name '*_R1_*' -exec rm '{}' \;
      find demultiplexed -name '*_R2_*' -exec rm '{}' \;
   fi
done
# Now, we trim assembled and unassembled reads separately as well. Hopefully,
# assembled reads should not need much trimming.
   for i in `seq 1 24`; do
      if [ ! -e trimmed/${SAMPLE[$i]}_R1.fastq ] && [ ! -e trimmed/${SAMPLE[$i]}_R1.fastq.gz ]; then
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

## Clean up.
for i in `seq 1 24`; do
   if [ ! -e demultiplexed/${SAMPLE[$i]}.fastq.gz ] && [ -e demultiplexed/${SAMPLE[$i]}.fastq ]; then
      gzip demultiplexed/${SAMPLE[$i]}.fastq &
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R1.fastq.gz ] && [ -e demultiplexed/${SAMPLE[$i]}_R1.fastq ]; then
      gzip demultiplexed/${SAMPLE[$i]}_R1.fastq &
   fi
   if [ ! -e demultiplexed/${SAMPLE[$i]}_R2.fastq.gz ] && [ -e demultiplexed/${SAMPLE[$i]}_R2.fastq ]; then
      gzip demultiplexed/${SAMPLE[$i]}_R2.fastq &
   fi
   if [ ! -e merged/${SAMPLE[$i]}.assembled.fastq.gz ] && [ -e merged/${SAMPLE[$i]}.assembled.fastq ]; then
      gzip merged/${SAMPLE[$i]}.assembled.fastq &
   fi
   if [ ! -e merged/${SAMPLE[$i]}.unassembled.forward.fastq.gz ] && [ -e merged/${SAMPLE[$i]}.unassembled.forward.fastq ]; then
      gzip merged/${SAMPLE[$i]}.unassembled.forward.fastq.gz &
   fi
   if [ ! -e merged/${SAMPLE[$i]}.unassembled.reverse.fastq.gz ] && [ -e merged/${SAMPLE[$i]}.unassembled.reverse.fastq ]; then
      gzip merged/${SAMPLE[$i]}.unassembled.reverse.fastq &
   fi
   if [ ! -e trimmed/${SAMPLE[$i]}.fastq.gz ] && [ -e trimmed/${SAMPLE[$i]}.fastq ]; then
      gzip trimmed/${SAMPLE[$i]}.fastq &
   fi
   if [ ! -e trimmed/${SAMPLE[$i]}_R1.fastq.gz ] && [ -e trimmed/${SAMPLE[$i]}_R1.fastq ]; then
      gzip trimmed/${SAMPLE[$i]}_R1.fastq &
   fi
   if [ ! -e trimmed/${SAMPLE[$i]}_R2.fastq.gz ] && [ -e trimmed/${SAMPLE[$i]}_R2.fastq ]; then
      gzip trimmed/${SAMPLE[$i]}_R2.fastq &
   fi
   wait
done

## Clean up.
#for i in `seq 1 24`; do
#   if [ -e merged/${SAMPLE[$i]}.assembled.fastq ]; then
#      rm demultiplexed/${SAMPLE[$i]}*.fastq
#      rm trimmed/${SAMPLE[$i]}*.fastq.gz
#   fi
#done
