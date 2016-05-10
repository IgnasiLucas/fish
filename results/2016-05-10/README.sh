#!/bin/bash
#
#				2016-05-10
#				----------
#
# I will start running again pyrad to cluster merged reads and build the
# consensus, using more threads here. I will monitor the progress and avort
# the process from 2016-05-03 if this one is expected to finish much earlier.
#

FASTQ=$(dirname $(pwd))"/2016-05-03"
if [ ! -e params.txt ]; then
   pyrad -n
                    #==** parameter inputs for pyRAD version 3.0.64  **======================== affected step ==
   sed -i '/## 6. /c\TAA,CGG                   ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)' params.txt
   sed -i '/## 7. /c\24                        ## 7. N processors (parallel)                           (all)'   params.txt
   sed -i '/## 8. /c\1                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)' params.txt
   sed -i '/## 9. /c\30                        ## 9. maxN: max number of Ns in reads                   (s2)'    params.txt
   sed -i '/## 10./c\.86                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)' params.txt
   sed -i '/## 11./c\ddrad                     ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(see docs)   (all)'   params.txt
                    #==== optional params below this line ===================================  affected step ==
   sed -i '/## 18./c\pooled.fastq.gz        ## 18.opt.: loc. of de-multiplexed data                 (s2)'       params.txt
   sed -i '/## 21./c\1                      ## 21.opt.: filter: def=0=NQual 1=NQual+adapt. 2=strict (s2)'       params.txt
   sed -i '/## 22./c\0.0000001,0.2          ## 22.opt.: a priori E,H (def=0.001,0.01, if not estim.)(s5)'       params.txt
   sed -i '/## 23./c\15                     ## 23.opt.: max N in consensus sequence (def. 5)        (s5)'       params.txt
   sed -i '/## 25./c\3                      ## 25.opt.: ploidy: max alleles in cons seq (def=2)     (s4,s5)'    params.txt
   sed -i '/## 31./c\6                      ## 31.opt.: maj. base call at depth <x (def.x=mindepth) (s5)'       params.txt
   sed -i '/## 32./c\50                     ## 32.opt.: keep trimmed reads (def=0). Enter min length(s2)'       params.txt
   sed -i '/## 37./c\32                     ## 37.opt.: vsearch max threads per job (def.=6)        (s3,s6)'    params.txt
fi

if [ ! -e pooled.fastq.gz ]; then
   ln -s $FASTQ/pooled.fastq.gz ./
fi

if [ ! -d edits ]; then mkdir edits; fi
if [ ! -e edits/pooled.edit ]; then
   ln -s $FASTQ/edits/pooled.edit edits/pooled.edit
fi
if [ ! -e edits/pooled.derep ]; then
   ln -s $FASTQ/edits/pooled.derep edits/pooled.derep
fi
if [ ! -d stats ]; then mkdir stats; fi
if [ ! -e stats/s2.rawedit.txt ]; then
   ln -s $FASTQ/stats/s2.rawedit.txt stats/s2.rawedit.txt
fi

touch checkpoints

if ! grep -q pyrad checkpoints; then
   echo pyrad >> checkpoints
   pyrad -p params.txt -s 35 &
fi

sleep 5m

H=0
while [ $H -le 12 ]; do
   gawk -v H=$H '{
      F[$1] = 1
      F[$2] = 1
   }END{
      for (f in F) N++
      print H "\t" N
   }' clust.86/pooled.u >> monitor.txt
   sleep 1h
   H=$(( $H + 1 ))
done
