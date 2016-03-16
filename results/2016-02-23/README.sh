#!/bin/bash
#
#				2016-02-23
#				----------
# After choosing a maximum number of Ns, I should check the effect of the clustering
# threshold on the number and size of within sample clusters. Aparently, the larger
# the number of Ns allowed in a sequence, the lower the clustering threshold should
# be, because Ns do not contribute similarity.

FASTQ=`pwd | sed "s_2016-02-23_2015-12-16b/demultiplexed_"`
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

# The strategy to generate several parameters files is the same as before: I make
# a generic one for the first value of the clustering threshold, and then modify it
# to get the rest.
if [ ! -e pe/clust74.txt ]; then
   cd pe
      pyrad -n
      sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
      sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
      sed -i '/## 9. /c\30                ## 9. maxN: max number of Ns in reads (s2)' params.txt
      sed -i '/## 10. /c\.74              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
      sed -i '/## 11. /c\pairddrad        ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
      sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
      sed -i '/## 14. /c\c74              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
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
      mv params.txt clust74.txt
   cd ..
fi

if [ ! -e se/clust74.txt ]; then
   cd se
      pyrad -n
      sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
      sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
      sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
      sed -i '/## 9. /c\30                ## 9. maxN: max number of Ns in reads (s2)' params.txt
      sed -i '/## 10. /c\.74              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
      sed -i '/## 11. /c\ddrad            ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
      sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
      sed -i '/## 14. /c\c74              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
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
      mv params.txt clust74.txt
   cd ..
fi

for i in 76 78 80 82 84 86 88 90 92 94 96 98 99; do
   sed "/## 10. /c\.$i                ## 10.Wclust: clustering threshold as a decimal (s3,s6)" pe/clust74.txt > pe/clust$i.txt
   sed -i "/## 14. /c\c$i                ## 14. Prefix name for final output (no spaces) (s7)" pe/clust$i.txt
   sed "/## 10. /c\.$i                ## 10.Wclust: clustering threshold as a decimal (s3,s6)" se/clust74.txt > se/clust$i.txt
   sed -i "/## 14. /c\c$i                ## 14. Prefix name for final output (no spaces) (s7)" se/clust$i.txt
done

cd pe
   if [ ! -e stats/s2.rawedit.txt ]; then
      pyrad -p clust74.txt -s 2 1> s2.log 2> s2.err
   fi
cd ../se
   if [ ! -e stats/s2.rawedit.txt ]; then
      pyrad -p clust74.txt -s 2 1> s2.log 2> s2.err
   fi
cd ..

for i in 74 76 78 80 82 84 86 88 90 92 94 96 98 99; do
   cd pe
      if [ ! -e stats/s3.clusters$i.txt ]; then
         pyrad -p clust$i.txt -s 3 1> clust$i.log 2> clust$i.err
         mv stats/s3.clusters.txt stats/s3.clusters$i.txt
      fi
   cd ../se
      if [ ! -e stats/s3.clusters$i.txt ]; then
         pyrad -p clust$i.txt -s 3 1> clust$i.log 2> clust$i.err
         mv stats/s3.clusters.txt stats/s3.clusters$i.txt
      fi
   cd ..
done

if [ -e pe/summary.txt ]; then
   gawk 'BEGIN{
      i = 1
   }(FNR == 1) {
      ClustTh = substr(FILENAME, 12, 2)
      THRESHOLD[i] = ClustTh
      i++
   }((FNR < 26) && (/^StC|^BlC/)){
      SAMPLES[$1] = 1
      ClustNum[ClustTh "," $1] = $2
      MeanDepth[ClustTh "," $1] = $3
      ClustNumD3[ClustTh "," $1] = $5
      MeanDepthD3[ClustTh "," $1] = $6
   }END{
      HEADER = "Sample"
      for (j = 1; j <= i; j++){
         HEADER = HEADER "\tClustNum_" THRESHOLD[j] "\tMeanDepth_" THRESHOLD[j] "\tClustNumD3_" THRESHOLD[j] "\tMeanDepthD3_" THRESHOLD[j]
      }
      print HEADER
   }' pe/stats/s3.clusters* > pe/summary.txt
fi
