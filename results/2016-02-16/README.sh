#!/bin/bash
#
#				2016-02-16
#				----------
#
# A promissing analysis is the quartet inference under the coalescent
# model, by Julia Chifman and Laura Kubatko. It is implemented in paup.
# I installed a demo version of paup with the SVDQuartets algorithm,
# and I run it successfully. Preliminar runs support the parallel speciation
# hypothesis with strong support. However, the data matrix I used is full
# of missing data, and I want to be more restrictive on the minimum number
# of samples required to use a locus. Previous settings, from 2016-01-17,
# requires only 4 samples. I need to re-run step 7 from pyrad, generate the
# nexus file, include the PAUP commands in the nexus and run it.

PYRADDIR=`pwd | sed 's_2016-02-16_2015-12-17/se_'`

# I need the clust.95 and stats directories for pyrad.
if [ ! -d clust.95 ]; then
   ln -s $PYRADDIR/clust.95 ./clust.95
fi

if [ ! -d stats ]; then mkdir stats; fi

# I will create new configuration files for pyrad. All them should
# request only the nexus output.
if [ ! -e params.txt ]; then
   cp $PYRADDIR/params.txt ./params.txt
   sed -i "/## 30./c\n               ## 30.opt.: Output formats... (s7)" params.txt
fi

# For the records, I will run pyrad with all mincov values between
# 2 and 20. Then, I will use the .loci outputs to count how many loci
# are covered at least so many times.
for i in `seq 2 20`; do
   j=`printf "%02u" $i`
   if [ ! -e params$j.txt ]; then
      sed    "/## 12. /c\\$i               ## 12. MinCov: min samples in a final locus (s7)" params.txt > params$j.txt
      sed -i "/## 14. /c\mincov$j          ## 14. Prefix name for final output (no spaces) (s7)" params$j.txt
   fi
done

for i in `seq 2 20`; do
   j=`printf "%02u" $i`
   touch stats/mincov$j.stats
   if [ ! -e outfiles/mincov$j.nex ]; then
      pyrad -p params$j.txt -s 7
      rm outfiles/mincov$j.excluded_loci
      rm outfiles/mincov$j.phy
   fi
done

# Below, I add the paup block only to nexus files that may be
# worth running. In order to remove the paup block, use the
# block:
for i in `seq 2 20`; do
   j=`printf "%02u" $i`
   sed -i -n '/begin PAUP/q; p' outfiles/mincov$j.nex
done

for i in `seq 4 20`; do
   j=`printf "%02u" $i`
   if ! grep -q "begin PAUP" outfiles/mincov$j.nex; then
      echo -e "\nbegin PAUP;"                          >> outfiles/mincov$j.nex
      echo "   taxpartition species=StCf: StCf0037-StCf0050, StCa: StCa0001-StCa0019, BlCa: BlCa0065-BlCa0108, BlCl: BlCl0091-BlCl0116;" >> outfiles/mincov$j.nex
      echo "   log file=mincov$j.log start stop;"      >> outfiles/mincov$j.nex
      echo "   SVDQuartets speciestree=yes partition=species evalQuartets=all bootstrap=no showScores=yes showSV=no keepQuartets=no;" >> outfiles/mincov$j.nex
      echo "   SVDQuartets speciestree=yes partition=species evalQuartets=all bootstrap=yes nreps=6000 nthreads=6; showScores=no" >> outfiles/mincov$j.nex
      echo "   describetrees 1 / plot=phylogram;"      >> outfiles/mincov$j.nex
      echo "   SaveTrees file=mincov$j.tree brlens;"   >> outfiles/mincov$j.nex
      echo "   log stop;"                              >> outfiles/mincov$j.nex
      echo "   quit;"                                  >> outfiles/mincov$j.nex
      echo -e "end;\n"                                 >> outfiles/mincov$j.nex
   fi
done

for i in `seq 4 20`; do
   j=`printf "%02u" $i`
   if [ ! -e outfiles/mincov$j.log ]; then
      paup outfiles/mincov$j.nex
   fi
done

if [ ! -e alltrees.nex ]; then
   echo "#NEXUS" > alltrees.nex
   echo "BEGIN TREES;"  >> alltrees.nex
   echo "   TRANSLATE"  >> alltrees.nex
   echo "      1 StCf," >> alltrees.nex
   echo "      2 BlCa," >> alltrees.nex
   echo "      3 StCa," >> alltrees.nex
   echo "      4 BlCl"  >> alltrees.nex
   echo "      ;"       >> alltrees.nex
   for i in `seq 10 20`; do
      j=`printf "%02u" $i`
      grep "^tree " outfiles/mincov$j.tree | gawk -v J=$j '{gsub(/PAUP_1/,"mincov" J); print}' >> alltrees.nex
   done
   echo "END;" >> alltrees.nex
fi

for i in `seq 4 20`; do
   j=`printf "%02u" $i`
   if [ ! -e outfiles/mincov$j.scores ]; then
      gawk -f scores.awk outfiles/mincov$j.log > outfiles/mincov$j.scores
   fi
done
