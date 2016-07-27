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
   sed -i '/## 8. /c\5                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)' params.txt
   sed -i '/## 9. /c\30                        ## 9. maxN: max number of Ns in reads                   (s2)'    params.txt
   sed -i '/## 10./c\.86                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)' params.txt
   sed -i '/## 11./c\ddrad                     ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(see docs)   (all)'   params.txt
   sed -i '/## 13./c30                         ## 13. MaxSH: max inds with shared hetero site          (s7)'    params.txt
                    #==== optional params below this line ===================================  affected step ==
   sed -i '/## 18./c\pooled.fastq.gz        ## 18.opt.: loc. of de-multiplexed data                 (s2)'       params.txt
   sed -i '/## 21./c\1                      ## 21.opt.: filter: def=0=NQual 1=NQual+adapt. 2=strict (s2)'       params.txt
   sed -i '/## 22./c\0.0000001,0.2          ## 22.opt.: a priori E,H (def=0.001,0.01, if not estim.)(s5)'       params.txt
   sed -i '/## 23./c\100                    ## 23.opt.: max N in consensus sequence (def. 5)        (s5)'       params.txt
   sed -i '/## 24./c\100                    ## 24.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)'  params.txt
   sed -i '/## 25./c\100                    ## 25.opt.: ploidy: max alleles in cons seq (def=2)     (s4,s5)'    params.txt
   sed -i '/## 31./c\1                      ## 31.opt.: maj. base call at depth <x (def.x=mindepth) (s5)'       params.txt
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

# After running successfully for over one month, I obtain this result in step 5:
#
# taxon        	nloci	f1loci	f2loci	nsites  	npoly	poly
# pooled	2601328	2601307	76475	33558428	159006	0.0047382
#
#    ## nloci = number of loci
#    ## f1loci = number of loci with >N depth coverage
#    ## f2loci = number of loci with >N depth and passed paralog filter
#    ## nsites = number of sites across f loci
#    ## npoly = number of polymorphic sites in nsites
#    ## poly = frequency of polymorphic sites
#
# Note that only 76475, out of 2601307 pass the paralog filter. This is apparently
# due to setting ploidy to 3. My intention was to prevent filtering of apparently
# triploid loci, because the pooled reads do contain several individuals, and
# more than 2 alleles are to be expected in some sites. It is clear that I underestimated
# the number of haplotypes that can be present in a sample of 24 diploid individuals.
# I should have used a ploidy of 48, at least. I will set it to 50 and re-run step
# 5 manually. Here, I just edit the script and params.txt file, as if this is what
# it acutally run.
#
# Finally, I realized that the mistake was to have a mindepth (parameter 7) lower than
# the lowcount (parameter 31). After setting mindepth to 6 and lowcount to 1, step 5
# run much faster (using the faster majority rule most of the time), and keeping most
# of the clusters.
