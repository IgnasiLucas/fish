#!/bin/bash
#
#				2016-01-13
#				----------
#
# Here, I run jmodeltest first to determine the optimum parameters of the
# molecular evolution models for the loci selected before, and I will run
# *BEAST. While jmodeltest can be run from the script, Beauti and BEAST
# may be run through their GUI.

PREVIOUS=`pwd | sed 's/2016-01-13/2016-01-11/'`
JMODELTEST=~/bin/jmodeltest-2.1.7/jModelTest.jar

if [ ! -e selected_loci_2.txt ]; then
   ln -s $PREVIOUS/selected_loci_2.txt ./selected_loci_2.txt
fi

if [ ! -e selected_loci_3.txt ]; then
   ln -s $PREVIOUS/selected_loci_3.txt ./selected_loci_3.txt
fi

if [ ! -d nexus ]; then mkdir nexus; fi

for i in `cut -f 1 selected_loci_2.txt`; do
   if [ ! -e $i ] && ! [[ $i =~ ^# ]]; then
      ln -s $PREVIOUS/$i $i
   fi
done

for i in `cut -f 1 selected_loci_3.txt`; do
   if [ ! -e $i ] && ! [[ $i =~ ^# ]]; then
      ln -s $PREVIOUS/$i $i
   fi
done

if [ ! -d jmodeltest ]; then mkdir jmodeltest; fi

for i in `ls -1 nexus`; do
   OUTFILE=jmodeltest/`basename -s .nex $i`.out
   if [ ! -e $OUTFILE ]; then
      java -jar $JMODELTEST -d nexus/$i -BIC -AIC -AICc -DT -f -o $OUTFILE
   fi
done

# BEAST produces the list of species trees, from which I can strip the branch lengths
# and annotations, to look at the topologies. BEAST is totally oriented to rooted trees.
# It would have been better to run MrBayes, which also implements multispecies coalescence.

if [ ! -e topologies.txt ]; then
   BURNIN=30
   SAMPLES=`wc -l species_*.trees | gawk 'END{print $1}'`
   gawk -v SAMPLES=$SAMPLES -v BURNIN=$BURNIN '((/^tree/) && (NR > SAMPLES * BURNIN / 100)){
      gsub(/:[^\),]+\)/,")",$4)
      gsub(/:[^\),]+,/,",",$4)
      gsub(/\[[^\]]+\]/,"",$4)
      F[$4]++
   }END{
      for (f in F) print f "\t" F[f]
   }' species_*.trees > topologies.txt
fi
