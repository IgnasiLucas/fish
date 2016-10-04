#!/bin/bash
#
#				2016-10-03
#				----------
#
# Here, I repeat the clustering of merged reads from all samples, to make an
# an assembly. I wanted to use ipyrad, supposed to run faster. However, it forces
# me to run step 4, the joint estimate of sequencing error and heterozygosity,
# which assumes diploidy. I resort to use pyrad again, like in 2016-05-10.
#
#
# -------------------------------------------------------------------------
# SETUP
# -------------------------------------------------------------------------

FASTQ=$(dirname $(pwd))"/2015-12-16b/joined"

SAMPLE=(ZeroNotUsed St0001 St0003        St0015 St0016 St0019 St0037 St0039
                    St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
                    Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

NEWNAME=(ZeroNotUsed StCa0001 StCa0003          StCa0015 StCa0016 StCa0019 StCf0037 StCf0039
                     StCf0043 StCf0044 StCf0049 StCf0050 BlCa0065 BlCa0076 BlCa0080 BlCa0083
                     BlCa0104 BlCa0108 BlCl0091 BlCl0093 BlCl0094 BlCl0095 BlCl0098 BlCl0116)

BARCODE=(ZeroNotUsed AGCTA TCGAG         ACCGAGT TGGGTGCC GTCTTGCG
         CAATCC GTTCAA AGCTA TCGAG TAACCTG ACCGAGT
         TGGGTGCC GTCTTGCG CAATCC GTTCAA AGCTA TCGAG
         TAACCTG ACCGAGT TGGGTGCC GTCTTGCG CAATCC GTTCAA)

touch checkpoints

# I need to uncompress the fastq files for sabre.
for i in `seq 1 23`; do
   for j in 1 2; do
      if ! grep -q ${NEWNAME[$i]}_R$j.fastq checkpoints; then
         echo ${NEWNAME[$i]}_R$j.fastq >> checkpoints
         cp $FASTQ/${SAMPLE[$i]}_R$j.fastq.gz ${NEWNAME[$i]}_R$j.fastq.gz
         gunzip ${NEWNAME[$i]}_R$j.fastq.gz
      fi
   done
done

# ---------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------

function demultiplex {
   if [ ! -d demultiplex ]; then
      echo "Please, make directory demultiplex"
      exit
   fi
   echo "$2 $1_demultiplexed_R1.fastq $1_demultiplexed_R2.fastq" > $1.barcode
   sabre pe -m 1 \
            -f $1_R1.fastq \
            -r $1_R2.fastq \
            -b $1.barcode \
            -u $1_unknown_R1.fastq \
            -w $1_unknown_R2.fastq
   echo $2 | gawk -v S=$1 '{print substr($1,2) " " S "_emultiplexed_R1.fastq " S "_emultiplexed_R2.fastq"}' > $1.arcode

   # Allow for some reads to miss the first base of the barcode.
   sabre pe -m 0 \
            -f $1_unknown_R1.fastq \
            -r $1_unknown_R2.fastq \
            -b $1.arcode \
            -u $1_uunknown_R1.fastq \
            -w $1_uunknown_R2.fastq

   if [ -s $1_emultiplexed_R1.fastq ]; then
      cat_seqs -o $1_dem_R1.fastq $1_demultiplexed_R1.fastq $1_emultiplexed_R1.fastq
      cat_seqs -o $1_dem_R2.fastq $1_demultiplexed_R2.fastq $1_emultiplexed_R2.fastq
      mv $1_dem_R1.fastq $1_demultiplexed_R1.fastq
      mv $1_dem_R2.fastq $1_demultiplexed_R2.fastq
   fi
   rm $1_emultiplexed_*
   rm $1_R1.fastq
   rm $1_R2.fastq
   rm $1_unknown*
   rm $1_uunknown*
   rm $1.barcode $1.arcode
}

function merge {
   if [ ! -d merged ]; then
      echo "Please make dir 'merged'"
      exit
   fi
   pear -f $1_demultiplexed_R1.fastq \
        -r $1_demultiplexed_R2.fastq \
        --output merged/$1 \
        --min-overlap 10 \
        --threads 16 \
        --memory 10G > merged/$1_pear.log
   # The following reverse-complementat of non-assembled second reads is done
   # to preserve the original orientation of second reads.
   vsearch --fastx_revcomp merged/$1.unassembled.reverse.fastq \
           --fastqout merged/$1.unassembled.reverse.reversed.fastq
   mv merged/$1.unassembled.reverse.reversed.fastq merged/$1.unassembled.reverse.fastq
}

function trim_pairs {
   BARCODEREV=`echo $2 | gawk -f revcomp.awk`
   cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A "TA"$BARCODEREV"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
            -o $1_trimmed_R1.fastq \
            -p $1_trimmed_R2.fastq \
            -m 50 \
            -q 5 \
            --pair-filter=both \
            merged/$1.unassembled.forward.fastq \
            merged/$1.unassembled.reverse.fastq
}

# ---------------------------------------------------------------------------
# DEMULTIPLEX
# ---------------------------------------------------------------------------

if [ ! -d demultiplex ]; then mkdir demultiplex; fi
for i in `seq 1 23`; do
   if ! grep -q ${NEWNAME[$i]}_demultiplexed checkpoints; then
      echo ${NEWNAME[$i]}_demultiplexed >> checkpoints
      demultiplex ${NEWNAME[$i]} ${BARCODE[$i]} > demultiplex/${NEWNAME[$i]}.log &
   fi
done
wait

if [ ! -e demultiplex/summary ]; then
   echo -e "Sample\tTotal\tFullMatch\tPartialMatch\tNoMatch" > demultiplex/summary
   gawk '(FNR == 2){
      TOTAL[substr(FILENAME,13,8)] = substr($5,2)
   }(FNR == 4){
      FULL[substr(FILENAME,13,8)] = substr($7,2)
   }(FNR == 13){
      PART[substr(FILENAME,13,8)] = substr($7,2)
   }(FNR == 15){
      NOPE[substr(FILENAME,13,8)] = substr($8,2)
   }END{
      for (f in TOTAL) {
         printf("%s\t%u\t%u (%.2f\%)\t%u (%.2f\%)\t%u (%.2f\%)\n", f, TOTAL[f], FULL[f], 100*FULL[f]/TOTAL[f], PART[f], 100*PART[f]/TOTAL[f], NOPE[f], 100*NOPE[f]/TOTAL[f])
      }
   }' demultiplex/*.log | sort -nk 2,2 >> demultiplex/summary
fi

# ---------------------------------------------------------------------------
# MERGE
# ---------------------------------------------------------------------------

if [ ! -d merged ]; then mkdir merged; fi
for i in `seq 1 23`; do
   if ! grep -q ${NEWNAME[$i]}_merged checkpoints; then
      echo ${NEWNAME[$i]}_merged >> checkpoints
      merge ${NEWNAME[$i]} > merged/${NEWNAME[$i]}.log &
   fi
   if [ -s merged/${SAMPLE[$i]}.assembled.fastq ]; then
      rm ${SAMPLE[$i]}_demultiplexed_R*.fastq
   fi
done
wait

# Now, there should be three files per sample in merge: *.assembled.fastq,
# *.unassembled.forward.fastq, and *.unassembled.reverse.fastq.

if [ ! -e merged/summary ]; then
   echo -e "#Sample\tTotal \tMerged  \tNotMerged\tDiscarded" > merged/summary
   gawk '(/^Assembled reads \./){
      gsub(/,/,""); A[substr(FILENAME,8,8)] = $4;  T[substr(FILENAME,8,8)] = $6; pA[substr(FILENAME,8,8)] = $7
   }(/^Discarded reads \./){
      gsub(/,/,""); D[substr(FILENAME,8,8)] = $4; pD[substr(FILENAME,8,8)] = $7
   }(/^Not assembled/){
      gsub(/,/,""); N[substr(FILENAME,8,8)] = $5; pN[substr(FILENAME,8,8)] = $8
   }END{
      for (f in T) print f "\t" T[f] "\t" A[f] " " pA[f] "\t" N[f] " " pN[f] "\t" D[f] " " pD[f]
   }' merged/*_pear.log | sort -k 1,1 >> merged/summary
fi

# ---------------------------------------------------------------------------
# TRIMMING
# ---------------------------------------------------------------------------
#
# I only trim paired ends.

if [ ! -d trimming ]; then mkdir trimming; fi

for i in `seq 1 23`; do
   if ! grep -q ${NEWNAME[$i]}_trimmed_paired checkpoints; then
      echo ${NEWNAME[$i]}_trimmed_paired >> checkpoints
      trim_pairs ${NEWNAME[$i]} ${BARCODE[$i]} > trimming/${NEWNAME[$i]}.log &
   fi
done
wait

if [ ! -e trimming/summary ]; then
   echo -e "#Sample\tPairs\tRead1Trimmed\tRead2Trimmed" > trimming/summary
   gawk '(/^Total read p/){
      gsub(/,/,""); T[substr(FILENAME,10,8)] = $5
   }(/Read 1 with/){
      gsub(/,/,""); A1[substr(FILENAME,10,8)] = $5 " " $6
   }(/Read 2 with/){
      gsub(/,/,""); A2[substr(FILENAME,10,8)] = $5 " " $6
   }END{
      for (f in T) print f "\t" T[f] "\t" A1[f] "\t" A2[f]
   }' trimming/*.log | sort -k 1,1 >> trimming/summary
fi

# ---------------------------------------------------------------------------
# PYRAD
# ---------------------------------------------------------------------------
#
# Here, I will repeat the clustering of assembled reads from all samples, following
# 2016-05-10.

if [ ! -e pooled.fastq.gz ]; then
   cat_seqs -o pooled.fastq.gz -z -t fastq merged/*.assembled.fastq
fi

if [ ! -e params.txt ]; then
   pyrad -n
                    #==** parameter inputs for pyRAD version 3.0.64  **======================== affected step ==
   sed -i '/## 6. /c\TAA,CGG                   ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)' params.txt
   sed -i '/## 7. /c\50                        ## 7. N processors (parallel)                           (all)'   params.txt
   sed -i '/## 8. /c\5                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)' params.txt
   sed -i '/## 9. /c\5                         ## 9. maxN: max number of Ns in reads                   (s2)'    params.txt
   sed -i '/## 10./c\.86                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)' params.txt
   sed -i '/## 11./c\ddrad                     ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(see docs)   (all)'   params.txt
   sed -i '/## 13./c30                         ## 13. MaxSH: max inds with shared hetero site          (s7)'    params.txt
                    #==== optional params below this line ===================================  affected step ==
   sed -i '/## 18./c\pooled.fastq.gz        ## 18.opt.: loc. of de-multiplexed data                 (s2)'       params.txt
   sed -i '/## 20./c\26                     ## 20.opt.: phred Qscore offset (def= 33)                    (s2)'  params.txt
   sed -i '/## 21./c\0                      ## 21.opt.: filter: def=0=NQual 1=NQual+adapt. 2=strict (s2)'       params.txt
   sed -i '/## 22./c\0.0000001,0.2          ## 22.opt.: a priori E,H (def=0.001,0.01, if not estim.)(s5)'       params.txt
   sed -i '/## 23./c\100                    ## 23.opt.: max N in consensus sequence (def. 5)        (s5)'       params.txt
   sed -i '/## 24./c\100                    ## 24.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)'  params.txt
   sed -i '/## 25./c\100                    ## 25.opt.: ploidy: max alleles in cons seq (def=2)     (s4,s5)'    params.txt
   sed -i '/## 31./c\3                      ## 31.opt.: maj. base call at depth <x (def.x=mindepth) (s5)'       params.txt
   sed -i '/## 32./c\50                     ## 32.opt.: keep trimmed reads (def=0). Enter min length(s2)'       params.txt
   sed -i '/## 37./c\50                     ## 37.opt.: vsearch max threads per job (def.=6)        (s3,s6)'    params.txt
fi

if ! grep -q pyrad checkpoints; then
   echo pyrad >> checkpoints
   pyrad -p params.txt -s 235
fi
