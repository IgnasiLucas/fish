#!/bin/bash
#
#				2016-02-22
#				----------
#
# One problem of using longer than usual reads is that they accumulate more
# low quality bases, which are masked as Ns. Also, the merging of the reads
# produces reads of different lengths, but all them are filtered by the number,
# not the proportion, of Ns. A higher number of Ns allowed should be accompanied
# by a lower clustering threshold.
#
# Here, I will check the effect of the maximum number of Ns in the number of
# reads that get filtered out.
#
# Note that I use here the demultiplexed, not trimmed, fastq files. I decided
# to override the trimming of adapter sequences done by cutadapt because: pyrad
# can also trim adapter sequences (probably not so sensitive, unless strict=2),
# and because the uneven trimming of members of a not merged pair results in
# problems with the not merged set.

FASTQ=`pwd | sed "s_2016-02-22_2015-12-16b/demultiplexed_"`
SAMPLE=(ZeroNotUsed St0001 St0003 St0006 St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)
NEWNAME=(ZeroNotUsed StCa0001 StCa0003 StCa0006 StCa0015 StCa0016 StCa0019 StCf0037 StCf0039
         StCf0043 StCf0044 StCf0049 StCf0050 BlCa0065 BlCa0076 BlCa0080 BlCa0083
         BlCa0104 BlCa0108 BlCl0091 BlCl0093 BlCl0094 BlCl0095 BlCl0098 BlCl0116)

if [ ! -d pe ]; then mkdir pe; fi
if [ ! -d pe/fastq ]; then mkdir pe/fastq; fi
if [ ! -d se ]; then mkdir se; fi
if [ ! -d se/fastq ]; then mkdir se/fastq; fi

# I need to specify not only the origin, (St or Bl) but either the species or
# the morphotype in the sample name.
for i in `seq 1 24`; do
   if [ ! -e pe/fastq/${NEWNAME[$i]}_R1.fastq ]; then
      ln -s $FASTQ/${SAMPLE[$i]}_R1.fastq pe/fastq/${NEWNAME[$i]}_R1.fastq
   fi
   if [ ! -e pe/fastq/${NEWNAME[$i]}_R2.fastq ]; then
      ln -s $FASTQ/${SAMPLE[$i]}_R2.fastq pe/fastq/${NEWNAME[$i]}_R2.fastq
   fi
   if [ ! -e se/fastq/${NEWNAME[$i]}.fastq ]; then
      ln -s $FASTQ/${SAMPLE[$i]}.fastq se/fastq/${NEWNAME[$i]}.fastq
   fi
done

if [ ! -e pe/maxN10.txt ]; then
   cd pe
      pyrad -n
      sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
      sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
      sed -i '/## 9. /c\10                ## 9. maxN: max number of Ns in reads (s2)' params.txt
      sed -i '/## 10. /c\.90              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
      sed -i '/## 11. /c\pairddrad        ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
      sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
      sed -i '/## 14. /c\N10              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
      sed -i '/## 17./c\St0006            ## 17.opt.: exclude taxa (list or prefix*) (s7)' params.txt
      sed -i '/## 18./c\fastq/*           ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
      sed -i '/## 21./c\1                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
      sed -i '/## 23./c\10                ## 23.opt.: max N in consensus sequence (def. 5) (s5)' params.txt
      sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
      sed -i '/## 32./c\75                ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)' params.txt
      sed -i '/## 36./c\64                ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
      echo "StCa 2 StCa*" >> params.txt
      echo "StCf 2 StCf*" >> params.txt
      echo "BlCa 2 BlCa*" >> params.txt
      echo "BlCl 2 BlCl*" >> params.txt
      mv params.txt maxN10.txt
   cd ..
fi

if [ ! -e se/maxN10.txt ]; then
   cd se
      pyrad -n
      sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
      sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
      sed -i '/## 9. /c\10                ## 9. maxN: max number of Ns in reads (s2)' params.txt
      sed -i '/## 10. /c\.90              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
      sed -i '/## 11. /c\ddrad            ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
      sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
      sed -i '/## 14. /c\N10              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
      sed -i '/## 17./c\St0006            ## 17.opt.: exclude taxa (list or prefix*) (s7)' params.txt
      sed -i '/## 18./c\fastq/*           ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
      sed -i '/## 21./c\1                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
      sed -i '/## 23./c\10                ## 23.opt.: max N in consensus sequence (def. 5) (s5)' params.txt
      sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
      sed -i '/## 32./c\75                ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)' params.txt
      sed -i '/## 36./c\64                ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
      echo "StCa 2 StCa*" >> params.txt
      echo "StCf 2 StCf*" >> params.txt
      echo "BlCa 2 BlCa*" >> params.txt
      echo "BlCl 2 BlCl*" >> params.txt
      mv params.txt maxN10.txt
   cd ..
fi

for i in 15 20 25 30 35 40 45 50 55 60 65 70; do
   sed "/## 9. /c\\$i                  ## 9. maxN: max number of Ns in reads (s2)" pe/maxN10.txt > pe/maxN$i.txt
   sed -i "/## 14. /c\N$i                 ## 14. Prefix name for final output (no spaces) (s7)" pe/maxN$i.txt
   sed -i "/## 23./c\\$i                ## 23.opt.: max N in consensus sequence (def. 5) (s5)" pe/maxN$i.txt
   sed "/## 9. /c\\$i                  ## 9. maxN: max number of Ns in reads (s2)" se/maxN10.txt > se/maxN$i.txt
   sed -i "/## 14. /c\N$i                 ## 14. Prefix name for final output (no spaces) (s7)" se/maxN$i.txt
   sed -i "/## 23./c\\$i                ## 23.opt.: max N in consensus sequence (def. 5) (s5)" se/maxN$i.txt
done

for i in 10 15 20 25 30 35 40 45 50 55 60 65 70; do
   cd pe
      if [ ! -e stats/s2.rawedit.maxN$i.txt ]; then
         pyrad -p maxN$i.txt -s 2 1> maxN$i.log 2> maxN$i.err
         mv stats/s2.rawedit.txt stats/s2.rawedit.maxN$i.txt
         rm edits/*
      fi
   cd ../se
      if [ ! -e stats/s2.rawedit.maxN$i.txt ]; then
         pyrad -p maxN$i.txt -s 2 1> maxN$i.log 2> maxN$i.err
         mv stats/s2.rawedit.txt stats/s2.rawedit.maxN$i.txt
         rm edits/*
      fi
   cd ..
done

if [ ! -e pe/summary.txt ]; then
   echo sample > pe/summary.txt
   head -n 25 pe/stats/s2.rawedit.maxN10.txt | tail -n 24 | sort -k 1,1 | cut -f 1 >> pe/summary.txt
   for i in 10 15 20 25 30 35 40 45 50 55 60 65 70; do
      echo maxN$i > pe/z1
      head -n 25 pe/stats/s2.rawedit.maxN$i.txt | \
      tail -n 24 | sort -k 1,1 | cut -f 5 >> pe/z1
      paste pe/summary.txt pe/z1 > pe/z2; mv pe/z2 pe/summary.txt
   done
   echo total > pe/z1
   head -n 25 pe/stats/s2.rawedit.maxN10.txt | tail -n 24 | sort -k 1,1 | cut -f 2 >> pe/z1
   paste pe/summary.txt pe/z1 > pe/z2; mv pe/z2 pe/summary.txt
   rm pe/z1
fi

if [ ! -e se/summary.txt ]; then
   echo sample > se/summary.txt
   head -n 25 se/stats/s2.rawedit.maxN10.txt | tail -n 24 | sort -k 1,1 | cut -f 1 >> se/summary.txt
   for i in 10 15 20 25 30 35 40 45 50 55 60 65 70; do
      echo maxN$i > se/z1
      head -n 25 se/stats/s2.rawedit.maxN$i.txt | \
      tail -n 24 | sort -k 1,1 | cut -f 5 >> se/z1
      paste se/summary.txt se/z1 > se/z2; mv se/z2 se/summary.txt
   done
   echo total > se/z1
   head -n 25 se/stats/s2.rawedit.maxN10.txt | tail -n 24 | sort -k 1,1 | cut -f 2 >> se/z1
   paste se/summary.txt se/z1 > se/z2; mv se/z2 se/summary.txt
   rm se/z1
fi

if [ ! -e maxN.png ]; then
   R --no-save < plot.R
fi
