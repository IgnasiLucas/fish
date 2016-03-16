#!/bin/bash
#
#				2015-12-17
#				----------
#
# Here, I will run pyrad separately in folders se and pe. Preliminar runs proved
# that pyrad does not handle well paired ends that have been trimmed, probably
# because sometimes one of the reads of a pair may be too short. Indeed, pyrad
# can detect and trim adapters as well. So, it should work with non-trimmed data.
# Below, I re-run pyrad from the demultiplexed fastq files, not the trimmed ones,
# both from merged and non-merged reads.

#FASTQ=`pwd | sed "s_2015-12-17_2015-12-16b/trimmed_"`
FASTQ=`pwd | sed "s_2015-12-17_2015-12-16b/demultiplexed_"`
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

if [ ! -e pe/params.txt ]; then
   # Paired ends are the reads that did not get assembled. Their trimming with
   # cutadapt must have been efficient, and I expect fewer Ns in these reads than
   # in single ends.
   cd pe
   pyrad -n
   sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
   sed -i '/## 7. /c\4                 ## 7. N processors (parallel) (all)' params.txt
   sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
   sed -i '/## 9. /c\10                ## 9. maxN: max number of Ns in reads (s2)' params.txt
   sed -i '/## 10. /c\.95              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
   sed -i '/## 11. /c\pairddrad        ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
   sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
   sed -i '/## 14. /c\c95              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
   sed -i '/## 17./c\St0006            ## 17.opt.: exclude taxa (list or prefix*) (s7)' params.txt
   sed -i '/## 18./c\fastq/*           ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
   sed -i '/## 21./c\1                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
   sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
   sed -i '/## 36./c\4                 ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
   echo "StCa 2 StCa*" >> params.txt
   echo "StCf 2 StCf*" >> params.txt
   echo "BlCa 2 BlCa*" >> params.txt
   echo "BlCl 2 BlCl*" >> params.txt
   cd ..
fi

if [ ! -e se/params.txt ]; then
   cd se
   pyrad -n
   sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
   sed -i '/## 7. /c\4                 ## 7. N processors (parallel) (all)' params.txt
   sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
   sed -i '/## 9. /c\30                ## 9. maxN: max number of Ns in reads (s2)' params.txt
   sed -i '/## 10. /c\.95              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
   sed -i '/## 11. /c\ddrad            ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
   sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
   sed -i '/## 14. /c\c95              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
   sed -i '/## 17./c\St0006            ## 17.opt.: exclude taxa (list or prefix*) (s7)' params.txt
   sed -i '/## 18./c\fastq/*           ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
   sed -i '/## 21./c\1                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
   sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
   sed -i '/## 36./c\4                 ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
   echo "StCa 2 StCa*" >> params.txt
   echo "StCf 2 StCf*" >> params.txt
   echo "BlCa 2 BlCa*" >> params.txt
   echo "BlCl 2 BlCl*" >> params.txt
   cd ..
fi

cd pe
   pyrad -p params.txt -s 234567 1> pyrad.log 2> pyrad.err &
cd ../se
   pyrad -p params.txt -s 234567 1> pyrad.log 2> pyrad.err &
cd ..

