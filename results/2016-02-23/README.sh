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

for i in pe se; do
   if [ -e $i/summary.txt ]; then
      gawk -f summary.awk $i/stats/s3.clusters* > $i/summary.txt
      R --no-save < plot.R
   fi
done

# In retrospect, I realize I should have used a subsample of the reads to study
# the effect of the clustering threshold on the number of clusters and mean depth.
# The relationship between clustering threshold and number of clusters is direct,
# as expected. The question is what the optimum threshold is. Since divergence is
# low, I expect the undetermined bases (Ns) to drive the split of clusters as the
# required similarity increases. Paralogs must exist at any similarity level that
# will spuriously be clustered together. The higher the required similarity level,
# the larger the number of paralogs correctly distinguished as different clusters.
# But this cannot be observed, and may have a smooth, minor effect.
#
# Non-merged (paired) reads are clustered using only the first read (see the pyrad
# manual). However, they were filtered for a maximum number of Ns concatenated. Since
# I expect most Ns of a pair to be in the second read, not used in clustering, the
# similarity among non-merged reads should be higher than that among merged reads,
# even though the length of first reads is shorter (300 bp) than the length of merged
# reads (~500).
#
# I need to know how the Ns are distributed.

for i in `seq 1 24`; do
   if [ ! -e pe/${NEWNAME[$i]}.Ns ]; then
      echo -e "# First\tBoth" > pe/${NEWNAME[$i]}.Ns
      gawk -f countNpe.awk pe/edits/${NEWNAME[$i]}.derep >> pe/${NEWNAME[$i]}.Ns
   fi

   if [ ! -e se/${NEWNAME[$i]}.Ns ]; then
      gawk -f countNse.awk se/edits/${NEWNAME[$i]}.derep > se/${NEWNAME[$i]}.Ns
   fi
done

if [ ! -e Ns.png ]; then
   R --no-save < plot02.R
fi

# The plot Ns.png shows that indeed the non-merged reads contain higher proportions
# of Ns. This was not expected if non-merged reads were just sequenced from long templates
# with the same average base quality than merged reads. I realize a large portion of
# non merged reads have not been merged because of the poor quality of the second or
# also the first read. This encourages me to try to cluster non-merged, first reads
# with merged reads.
#
# Non-merged, first reads can have up to 10% of Ns, while merged reads rarely go
# above 5%. This percentages coincide roughly with the (di)similarity thresholds at
# which the number of clusters start growing: clustering thresholds above 90% similarity
# make the number of clusters grow much faster among non-merged reads; while the
# fast increase in the number of clusters of merged reads happens only around the
# 95% similarity threshold. This confirms the Ns (not paralogs) are driving the split
# of clusters with increasing similarity thresholds.
#
# This amounts to argue that it is safe to use a relatively low similarity threshold,
# because I don't expect a large number of clusters to be affected by paralogs.
#
# Sibelle Vilaça called my attention about missing reads: while step 3 in pyrad is
# not a filtering step, counting the number of reads in .derep and .clustS files
# shows a number of reads are missing. Curiously, the lower the similarity threshold,
# the larger the number of reads missing. This explains why the mean depth of clusters
# peak around thresholds 86% (non-merged) or 92% (merged), right before dropping. Those
# peaks were not expected if the number of reads being clustered were always the same.
#
# The fact is that step 3 filters out some reads that are first assigned to a cluster,
# and then discarded on the grounds of having too many gaps. The number of gaps in the
# pairwise comparison performed by vsearch is reported in .u files. The distributions
# may give hints on the optimum value of the maximum number of gaps allowed.
