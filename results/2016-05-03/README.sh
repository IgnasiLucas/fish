#!/bin/bash
#
#				2016-05-03
#				----------
#
# Here, I will demultiplex, merge, and trim the reads again, starting from the
# reads joined from the two sequencing runs. Then, I will create a reference
# with the consensus sequences from the merged reads and map the both merged
# and not merged reads to those consensus. The single reads that fail to map will
# be processed separately with pyrad.
#
# Since the analysis on 2015-12-16b, I've learned that both bbmerge.sh (from bbmap
# package) and vsearch can merge overlapping read pairs. But I decide to keep using
# PEAR, which has proved to be quite sensitive. For trimming, I also have some options,
# but I also stick to cutadapt, which gives a nice output and works just fine.
#
# For demultiplexing, I choose Sabre, over pyrad (stacks is another option) because
# it is easier, despite of the fact that I like pyrad's output better. Remember that
# this dataset had already been demultiplexed using Illumina indices. The current
# demultiplexing step is meant to filter reads with a wrong barcode, potentially
# coming from chimeric fragments.
#
# I exclude sample St0006 from the analysis.
#
# Below, I use a checkpoints-based conditional execution, to make it easier to
# delete files after being used.
#

# -------------------------------------------------------------------------
# SETUP
# -------------------------------------------------------------------------

FASTQ=$(dirname $(pwd))"/2015-12-16b/joined"
SAMPLE=(ZeroNotUsed St0001 St0003        St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)
BARCODE=(ZeroNotUsed AGCTA TCGAG         ACCGAGT TGGGTGCC GTCTTGCG
         CAATCC GTTCAA AGCTA TCGAG TAACCTG ACCGAGT
         TGGGTGCC GTCTTGCG CAATCC GTTCAA AGCTA TCGAG
         TAACCTG ACCGAGT TGGGTGCC GTCTTGCG CAATCC GTTCAA)

touch checkpoints

# I need to uncompress the fastq files for sabre.
for i in `seq 1 23`; do
   for j in 1 2; do
      if ! grep -q ${SAMPLE[$i]}_R$j.fastq checkpoints; then
         echo ${SAMPLE[$i]}_R$j.fastq >> checkpoints
         cp $FASTQ/${SAMPLE[$i]}_R$j.fastq.gz ${SAMPLE[$i]}_R$j.fastq.gz
         gunzip ${SAMPLE[$i]}_R$j.fastq.gz
      fi
   done
done

# ---------------------------------------------------------------------------
# DEMULTIPLEX
# ---------------------------------------------------------------------------

for i in `seq 1 23`; do
   if ! grep -q ${SAMPLE[$i]}_demultiplexed checkpoints; then
      echo ${SAMPLE[$i]}_demultiplexed >> checkpoints
      echo "${BARCODE[$i]} ${SAMPLE[$i]}_demultiplexed_R1.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq" > ${SAMPLE[$i]}.barcode
      sabre pe -m 1 \
               -f ${SAMPLE[$i]}_R1.fastq \
               -r ${SAMPLE[$i]}_R2.fastq \
               -b ${SAMPLE[$i]}.barcode \
               -u ${SAMPLE[$i]}_unknown_R1.fastq \
               -w ${SAMPLE[$i]}_unknown_R2.fastq > ${SAMPLE[$i]}_sabre.log
      echo ${BARCODE[$i]} | gawk -v S=${SAMPLE[$i]} '{print substr($1,2) " " S "_emultiplexed_R1.fastq " S "_emultiplexed_R2.fastq"}' > ${SAMPLE[$i]}.arcode

      # Allow for some reads to miss the first base of the barcode.
      sabre pe -m 0 \
               -f ${SAMPLE[$i]}_unknown_R1.fastq \
               -r ${SAMPLE[$i]}_unknown_R2.fastq \
               -b ${SAMPLE[$i]}.arcode \
               -u ${SAMPLE[$i]}_uunknown_R1.fastq \
               -w ${SAMPLE[$i]}_uunknown_R2.fastq > ${SAMPLE[$i]}_ssabre.log
      if [ -s ${SAMPLE[$i]}_emultiplexed_R1.fastq ]; then
         cat_seqs -o ${SAMPLE[$i]}_dem_R1.fastq ${SAMPLE[$i]}_demultiplexed_R1.fastq ${SAMPLE[$i]}_emultiplexed_R1.fastq
         cat_seqs -o ${SAMPLE[$i]}_dem_R2.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq ${SAMPLE[$i]}_emultiplexed_R2.fastq
         mv ${SAMPLE[$i]}_dem_R1.fastq ${SAMPLE[$i]}_demultiplexed_R1.fastq
         mv ${SAMPLE[$i]}_dem_R2.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq
      fi
      rm ${SAMPLE[$i]}_emultiplexed_*
      rm ${SAMPLE[$i]}_R1.fastq
      rm ${SAMPLE[$i]}_R2.fastq
      rm ${SAMPLE[$i]}_unknown*
      rm ${SAMPLE[$i]}_uunknown*
      rm ${SAMPLE[$i]}.barcode ${SAMPLE[$i]}.arcode
      if [ ! -d demultiplex ]; then mkdir demultiplex; fi
      mv *sabre.log demultiplex/
   fi
done

if [ ! -e demultiplex/summary ]; then
   gawk '(/^Total/){
      T[substr(FILENAME,13,6)] = substr($5,2)
   }(/for barcode/){
      M[substr(FILENAME,13,6)] = substr($7,2)
   }(/no barcode/){
      N[substr(FILENAME,13,6)] = substr($8,2)
   }END{
      for (f in T) print f "\t" T[f] + 0 "\t" M[f] + 0 "\t" N[f] + 0
   }' demultiplex/*_sabre.log | sort -k 1,1 > z1

   gawk '(/^Total/){
      T[substr(FILENAME,13,6)] = substr($5,2)
   }(/for barcode/){
      M[substr(FILENAME,13,6)] = substr($7,2)
   }(/no barcode/){
      N[substr(FILENAME,13,6)] = substr($8,2)
   }END{
      for (f in T) print f "\t" T[f] + 0 "\t" M[f] + 0 "\t" N[f] + 0
   }' demultiplex/*_ssabre.log | sort -k 1,1 > z2

   echo -e "#Sample\tTotal\tFullMatch\tPartialMatch\tNoMatch" | tee demultiplex/summary
   paste z1 z2 | gawk '{
      S = $1; T = $2; F = $3; pF = 100 * F / T; P = $7; pP = 100 * P / T; N = $8; pN = 100 * N / T
      printf("%s\t%u\t%u (%.2f\%)\t%u (%.2f\%)\t%u (%.2f\%)\n", S, T, F, pF, P, pP, N, pN)
   }' | tee -a demultiplex/summary
   rm z1 z2
fi

# --------------------------------------------------------------------------
# MERGE
# --------------------------------------------------------------------------

if [ ! -d merged ]; then mkdir merged; fi
for i in `seq 1 23`; do
   if ! grep -q ${SAMPLE[$i]}_merged checkpoints; then
      echo ${SAMPLE[$i]}_merged >> checkpoints
      pear -f ${SAMPLE[$i]}_demultiplexed_R1.fastq \
           -r ${SAMPLE[$i]}_demultiplexed_R2.fastq \
           --output merged/${SAMPLE[$i]} \
           --min-overlap 10 \
           --threads 16 \
           --memory 10G > merged/${SAMPLE[$i]}_pear.log
      # The following reverse-complementat of non-assembled second reads is done
      # to preserve the original orientation of second reads.
      vsearch --fastx_revcomp merged/${SAMPLE[$i]}.unassembled.reverse.fastq \
              --fastqout merged/${SAMPLE[$i]}.unassembled.reverse.reversed.fastq
      mv merged/${SAMPLE[$i]}.unassembled.reverse.reversed.fastq merged/${SAMPLE[$i]}.unassembled.reverse.fastq

      if [ -s merged/${SAMPLE[$i]}.assembled.fastq ]; then
         rm ${SAMPLE[$i]}_demultiplexed_R*.fastq
      fi
   fi
done

# Now, there should be three files per sample in merge: *.assembled.fastq,
# *.unassembled.forward.fastq, and *.unassembled.reverse.fastq.

if [ ! -e merged/summary ]; then
   echo -e "#Sample\tTotal \tMerged  \tNotMerged\tDiscarded" | tee merged/summary
   gawk '(/^Assembled reads \./){
      gsub(/,/,""); A[substr(FILENAME,8,6)] = $4;  T[substr(FILENAME,8,6)] = $6; pA[substr(FILENAME,8,6)] = $7
   }(/^Discarded reads \./){
      gsub(/,/,""); D[substr(FILENAME,8,6)] = $4; pD[substr(FILENAME,8,6)] = $7
   }(/^Not assembled/){
      gsub(/,/,""); N[substr(FILENAME,8,6)] = $5; pN[substr(FILENAME,8,6)] = $8
   }END{
      for (f in T) print f "\t" T[f] "\t" A[f] " " pA[f] "\t" N[f] " " pN[f] "\t" D[f] " " pD[f]
   }' merged/*_pear.log | sort -k 1,1 | tee -a merged/summary
fi


# ----------------------------------------------------------------------------
# TRIM
# ----------------------------------------------------------------------------

if [ ! -d trimmed ]; then mkdir trimmed; fi
for i in `seq 1 23`; do
   if ! grep -q ${SAMPLE[$i]}_trimmed_paired checkpoints; then
      echo ${SAMPLE[$i]}_trimmed_paired >> checkpoints
      BARCODEREV=`echo ${BARCODE[$i]} | gawk -f revcomp.awk`
      cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               -A "TA"$BARCODEREV"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
               -o ${SAMPLE[$i]}_trimmed_R1.fastq \
               -p ${SAMPLE[$i]}_trimmed_R2.fastq \
               -m 50 \
               -q 5 \
               --pair-filter=both \
               merged/${SAMPLE[$i]}.unassembled.forward.fastq \
               merged/${SAMPLE[$i]}.unassembled.reverse.fastq > trimmed/${SAMPLE[$i]}_cutadapt.pe.log
      if [ -s ${SAMPLE[$i]}_trimmed_R1.fastq ]; then
         rm merged/${SAMPLE[$i]}.unassembled.forward.fastq
         rm merged/${SAMPLE[$i]}.unassembled.reverse.fastq
      fi
   fi

   if ! grep -q ${SAMPLE[$i]}_trimmed_merged checkpoints; then
      echo ${SAMPLE[$i]}_trimmed_merged >> checkpoints
      cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT${BARCODE[$i]} \
               -o ${SAMPLE[$i]}_trimmed.fastq \
               -m 50 \
               merged/${SAMPLE[$i]}.assembled.fastq > trimmed/${SAMPLE[$i]}_cutadapt.se.log
      if [ -s ${SAMPLE[$i]}_trimmed.fastq ]; then
         rm merged/${SAMPLE[$i]}.assembled.fastq
      fi
   fi
done

if [ ! -e trimmed/summary ]; then
   echo -e "#Sample\tMerged\tMergedTrimmed\tMergedShort" > z1
   gawk '(/^Total reads/){
      gsub(/,/,""); T[substr(FILENAME,9,6)] = $4
   }(/^Reads with adap/){
      gsub(/,/,""); A[substr(FILENAME,9,6)] = $4 " " $5
   }(/^Reads that were/){
      gsub(/,/,""); S[substr(FILENAME,9,6)] = $6 " " $7
   }END{
      for (f in T) print f "\t" T[f] "\t" A[f] "\t" S[f]
   }' trimmed/*.se.log | sort -k 1,1 >> z1

   echo -e "#Pairs\tRead1Trimmed\tRead2Trimmed" > z2
   gawk '(/^Total read p/){
      gsub(/,/,""); T[substr(FILENAME,9,6)] = $5
   }(/Read 1 with/){
      gsub(/,/,""); A1[substr(FILENAME,9,6)] = $5 " " $6
   }(/Read 2 with/){
      gsub(/,/,""); A2[substr(FILENAME,9,6)] = $5 " " $6
   }END{
      for (f in T) print f "\t" T[f] "\t" A1[f] "\t" A2[f]
   }' trimmed/*.pe.log | sort -k 1,1 >> z2
   cut -f 2,3,4 z2 > z3
   paste z1 z3 | tee trimmed/summary
   rm z1 z2 z3
fi

# ----------------------------------------------------------
# DEREPLICATE & CLUSTER
# ----------------------------------------------------------
#
# While I could use vsearch directly to cluster all merged reads
# and obtain their consensus sequences, I notice that pyrad spends
# a large number of lines of code to replace bases with phred qualities
# lower than 20 into Ns. Then, de-replication must reduce a lot the
# number of reads, and speed up clustering. Let's use pyrad.
#
# I pool all merged fastq files, as if they were a single file, so that
# a single consensus sequence is generated for each locus on step 5. The
# problem is that a pooled dataset is not diploid: more than two alleles
# could appear in a site (parameter 25), virtually every variable site
# in a locus will look like heterozygous (parameter 24), and the number
# of SNPs might be high (parameter 26). On top of that, the estimation
# of heterozygosity and error rate, which are very time-consuming, would
# be biased: low MAF sites would look like errors, and only high-MAF
# sites would look like heterozygous. I could set the error rate and the
# heterozygosity to favor a variable consensus, with a low E and a high H.
# However, the binomial distribution of reads of each allele expected in
# a heterozygous site would make most sites with low MAF look like errors.

if ! grep -q dereplicate checkpoints; then
   echo dereplicate >> checkpoints
   if [ ! -e pooled.fastq.gz ]; then
      cat_seqs -o pooled.fastq.gz -z -t fastq *_trimmed.fastq
   fi
   pyrad -n
                    #==** parameter inputs for pyRAD version 3.0.64  **======================== affected step ==
   sed -i '/## 6. /c\TAA,CGG                   ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)' params.txt
   sed -i '/## 7. /c\24                        ## 7. N processors (parallel)                           (all)'   params.txt
   sed -i '/## 8. /c\1                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)' params.txt
   sed -i '/## 9. /c\30                        ## 9. maxN: max number of Ns in reads                   (s2)'    params.txt
   sed -i '/## 10./c\.86                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)' params.txt
   sed -i '/## 11./c\ddrad                     ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(see docs)   (all)'   params.txt
   sed -i '/## 18./c\pooled.fastq.gz           ## 18.opt.: loc. of de-multiplexed data                 (s2)'    params.txt
   sed -i '/## 21./c\1                         ## 21.opt.: filter: def=0=NQual 1=NQual+adapt. 2=strict (s2)'    params.txt
   sed -i '/## 22./c\0.0000001,0.2             ## 22.opt.: a priori E,H (def=0.001,0.01, if not estim.)(s5)'    params.txt
   sed -i '/## 23./c\15                        ## 23.opt.: max N in consensus sequence (def. 5)        (s5)'    params.txt
   sed -i '/## 25./c\3                         ## 25.opt.: ploidy: max alleles in cons seq (def=2)     (s4,s5)' params.txt
   sed -i '/## 31./c\6                         ## 31.opt.: maj. base call at depth <x (def.x=mindepth) (s5)'    params.txt
   sed -i '/## 32./c\50                        ## 32.opt.: keep trimmed reads (def=0). Enter min length(s2)'    params.txt
   sed -i '/## 37./c\24                        ## 37.opt.: vsearch max threads per job (def.=6)        (s3,s6)' params.txt

   pyrad -p params.txt -s 235
fi

# After several days running, a plot of the number of different reads
# included in the clust.86/pooled.u file is shown to be a perfectly linear
# function of time, and quite slow: 3194 new reads per hour. At this pace,
# it would take more than 90 days to finish. 
