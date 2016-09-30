#!/bin/bash
#
#				2016-09-19
#				----------
# In 2016-02-23, I checked the effect of the clustering threshold on the number and size
# the clusters produced, holding the maximum number of Ns per read constant at 30. After
# the last results from 2016-09-14, I suspect that the clustering from 2016-05-10 fails
# to cluster all the reads from the same locus. Splitting loci in separate clusters is
# unexpected with a similarity threshold as low as 86%. Unless the presence of Ns makes
# the threshold more stringent than it seams. The solution cannot be to keep lowering
# the similarity threshold, but to limit the number of Ns allowed.
#
# Here, I want to check the effect of the maximum number of Ns on the number and size
# of the clusters. I will use only 1 sample, since they all behave similarly. I choose
# BlCl0091.


FASTQ=$(dirname $(pwd))"/2016-05-03"
if [ ! -d pe ]; then mkdir pe; fi
if [ ! -d pe/fastq ]; then mkdir pe/fastq; fi
if [ ! -e pe/fastq/BlCl0091_R1.fastq ]; then
   ln -s $FASTQ/Bl0091_trimmed_R1.fastq pe/fastq/BlCl0091_R1.fastq
fi
if [ ! -e pe/fastq/BlCl0091_R2.fastq ]; then
   ln -s $FASTQ/Bl0091_trimmed_R2.fastq pe/fastq/BlCl0091_R2.fastq
fi
if [ ! -d se ]; then mkdir se; fi
if [ ! -d se/fastq ]; then mkdir se/fastq; fi
if [ ! -e se/fastq/BlCl0091.fastq ]; then
   ln -s $FASTQ/Bl0091_trimmed.fastq se/fastq/BlCl0091.fastq
fi

# pe/
#    fastq/
#    N25/
#        params.txt
#        stats
#        edits
#    N20/
#        ...
#    N15/
#        ...
#    ...
# se/
#    fastq/
#    N25/
#        ...
#    ...

if [ ! -d pe/N25 ]; then mkdir pe/N25; fi

if [ ! -e pe/N25/params.txt ]; then
   pyrad -n
   sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
   sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
   sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
   sed -i '/## 9. /c\25                ## 9. maxN: max number of Ns in reads (s2)' params.txt
   sed -i '/## 10. /c\.86              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
   sed -i '/## 11. /c\pairddrad        ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
   sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
   sed -i '/## 14. /c\cN25             ## 14. Prefix name for final output (no spaces) (s7)' params.txt
   sed -i '/## 18./c\../fastq/*        ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
   sed -i '/## 21./c\0                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
   sed -i '/## 23./c\10                ## 23.opt.: max N in consensus sequence (def. 5) (s5)' params.txt
   sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
   sed -i '/## 32./c\75                ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)' params.txt
   sed -i '/## 36./c\6                 ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
   mv params.txt pe/N25/params.txt
fi

if [ ! -d se/N25 ]; then mkdir se/N25; fi
if [ ! -e se/N25/params.txt ]; then
   pyrad -n
   sed -i '/## 6. /c\TAA,CGG           ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG) (s1,s2)' params.txt
   sed -i '/## 7. /c\24                ## 7. N processors (parallel) (all)' params.txt
   sed -i '/## 8. /c\4                 ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
   sed -i '/## 9. /c\25                ## 9. maxN: max number of Ns in reads (s2)' params.txt
   sed -i '/## 10. /c\.86              ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
   sed -i '/## 11. /c\ddrad            ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
   sed -i '/## 13. /c\6                ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
   sed -i '/## 14. /c\N25              ## 14. Prefix name for final output (no spaces) (s7)' params.txt
   sed -i '/## 18./c\../fastq/*        ## 18.opt.: loc. of de-multiplexed data (s2)' params.txt
   sed -i '/## 21./c\0                 ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
   sed -i '/## 23./c\10                ## 23.opt.: max N in consensus sequence (def. 5) (s5)' params.txt
   sed -i '/## 30./c\*                 ## 30.opt.: Output formats... (s7)' params.txt
   sed -i '/## 32./c\75                ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)' params.txt
   sed -i '/## 36./c\6                 ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
   mv params.txt se/N25/params.txt
fi

for i in 20 15 10 5; do
   if [ ! -d pe/N$i ]; then mkdir pe/N$i; fi
   sed "/## 9. /c\\$i                  ## 9.maxN: max number of Ns in reads (s2)" pe/N25/params.txt > pe/N$i/params.txt
   sed -i "/## 14. /c\N$i                ## 14. Prefix name for final output (no spaces) (s7)" pe/N$i/params.txt
   if [ ! -d se/N$i ]; then mkdir se/N$i; fi
   sed "/## 9. /c\\$i                  ## 9.maxN: max number of Ns in reads (s2)" se/N25/params.txt > se/N$i/params.txt
   sed -i "/## 14. /c\N$i                ## 14. Prefix name for final output (no spaces) (s7)" se/N$i/params.txt
done

for i in 25 20 15 10 5; do
   cd pe/N$i
      if [ ! -e stats/s3.clusters.txt ]; then
         pyrad -p params.txt -s 23 1> pyrad.log 2> pyrad.err &
      fi
   cd ../../se/N$i
      if [ ! -e stats/s3.clusters.txt ]; then
         pyrad -p params.txt -s 23 1> pyrad.log 2> pyrad.err &
      fi
   cd ../..
done
wait

if [ ! -e summary.txt ]; then
   echo -e "Set\tMaxN\tReads\tPassed\tClusters\tDepth" > summary.txt
   find ./pe -name 's2.rawedit.txt' -exec gawk '(/^BlCl0091/){
      split(FILENAME,A,"/")
      print "pe\t" substr(A[3],2) "\t" $2 "\t" $5
   }' '{}' \; | sort -nk 2 > z1

   find ./se -name 's2.rawedit.txt' -exec gawk '(/^BlCl0091/){
      split(FILENAME,A,"/")
      print "se\t" substr(A[3],2) "\t" $2 "\t" $5
   }' '{}' \; | sort -nk 2 >> z1

   find ./pe -name 's3.clusters.txt' -exec gawk '(/^BlCl0091/){
      split(FILENAME,A,"/")
      print substr(A[3],2) "\t" $2 "\t" $3
   }' '{}' \; | sort -nk 1 | cut -f 2,3 > z2

   find ./se -name 's3.clusters.txt' -exec gawk '(/^BlCl0091/){
      split(FILENAME,A,"/")
      print substr(A[3],2) "\t" $2 "\t" $3
   }' '{}' \; | sort -nk 1 | cut -f 2,3 >> z2

   paste z1 z2 >> summary.txt
   rm z1 z2
fi

if [ ! -e reads.png ] || [ ! -e clusters.png ] || [ ! -e depth.png ]; then
   R --save < plot.R
fi

# CONCLUSIONS
# -----------
#
# The lack of saturation of the number of clusters as the number of reads
# increases (depth.png), or the independence of depth with respect to the
# number of clusters, strongly suggests that there are just not enough
# sequences. Not all loci got sequenced in each sample. However, something
# else could contribute to this pattern. It is clear that many reads are
# of low quality and contain several undetermined bases (reads.png). Given
# that previous runs of the clustering step produced large numbers of
# highly related clusters (see 2016-09-14), it is possible that the abundance
# of Ns are making many clusters be split spuriously.

