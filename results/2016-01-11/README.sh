#!/bin/bash
#
#				2016-01-11
#				----------
#
# Up to now, only the pyrad analysis of the assembled reads produced some useful
# data. Before optimizing the use of paired ends, I need to set up the next step
# of the analysis. Some preliminar runs of *BEAST makes it clear that the loci
# sequenced in at least one individual of the four populations contain noise. I
# suspect some variable sites are indeed sequencing errors. This is a big problem
# because for the coalescence analysis I will be selecting sites with more than
# one variable (phylogenetically informative) site.
#
# In order to minimize the noise, I plan to check the distribution of variable
# sites along the assembled reads and the type of variation (transition or trans-
# version).
#
# Assembled reads are mostly made of two 300 bp reads. There are 3 types of regions
# in an assembly: covered by the 1st read (higher quality than the second), covered
# by both reads (overlap), and covered by the second read only (lower quality).
# The quality of the overlap should be high. The length of each part is variable.
# Although not all reads are 300 bp, most are. So, if assembled reads were shorter
# than 300 bp, it is probably because the sequenced fragment was shorter as well,
# and the whole assembly must be covered by both reads, with excellent quality.
# On the other hand, assemblies longer than 300 bp must consist on L-300 bp covered
# only by read 1, 600-L bp of overlap, and other L-300 bp covered only by read 2.
#
# The script VariationDistribution.awk will count the number of variable sites per
# region, as well as the number of sites and sequences in each region.

PREVIOUS=`pwd | sed 's/2016-01-11/2015-12-17/'`

if [ ! -e c95.loci ]; then
   ln -s $PREVIOUS/se/outfiles/c95.loci c95.loci
fi

if [ ! -e c95.alleles ]; then
   ln -s $PREVIOUS/se/outfiles/c95.alleles c95.alleles
fi

if [ ! -e variantPositions.txt ]; then
   gawk -f VariationDistribution.awk c95.loci > variantPositions.txt
fi

# The results do not show any obvious excess of variant sites in regions covered
# only by read 2, which is of lower quality. A simpler way to look at it would be
# to check the distribution of relative positions of variable sites.

if [ ! -e lengths.txt ]; then
   gawk '(/^\/\//){print LENGTH}(/^>/){LENGTH = length($2)}' c95.loci > lengths.txt
fi

if [ ! -e lengths.png ]; then
   R --no-save < plot1.R
fi

if [ ! -e relativePositions_0-200.txt ]; then
   gawk -v MAX=200 -f RelativePositions.awk c95.loci > relativePositions_0-200.txt
fi

if [ ! -e relativePositions_201-300.txt ]; then
   gawk -v MIN=201 -v MAX=300 -f RelativePositions.awk c95.loci > relativePositions_201-300.txt
fi

if [ ! -e relativePositions_301-400.txt ]; then
   gawk -v MIN=301 -v MAX=400 -f RelativePositions.awk c95.loci > relativePositions_301-400.txt
fi

if [ ! -e relativePositions_401-500.txt ]; then
   gawk -v MIN=401 -v MAX=500 -f RelativePositions.awk c95.loci > relativePositions_401-500.txt
fi

if [ ! -e relativePositions_501-600.txt ]; then
   gawk -v MIN=501 -v MAX=600 -f RelativePositions.awk c95.loci > relativePositions_501-600.txt
fi

if [ ! -e relativePositions.png ]; then
   R --no-save < plot2.R
fi

# The bulk of the loci are more than 500 bp long. Those accumulate an excess of
# variable sites right after the middle, between the midpoint and the 60% of the
# length. This has a straighforward interpretation: they are the last positions
# of read 2 before the overlap region. More sequencing errors are expected there.
# Thus, it makes sense to remove those positions from the analysis. I can substitute
# them by Ns.
#
# Shorter sequences are not so easy to understand or correct. Below 300 bp, the last
# 0.5% of the assembly accumulates the errors and could be trimmed. However, the reason
# of this bias is unclear, and the absolute number of loci (let alone informative sites)
# recovered, relatively low.
#
# For the moment, I'll focus on the longer, more promissing loci.

if [ ! -d nexus ]; then
   mkdir nexus
   gawk -v MINLENGTH=500 -v STARTRM=0.5 -v ENDRM=0.58 -f makenexus_trim.awk c95.alleles
   mv *.nex nexus/
fi

# Below, I select loci with only 2 informative sites and write an haplotypes file
# with the haplotypes in each locus.

if [ ! -e haplotypes ]; then
   gawk '((FILENAME ~ /_I2_/) && (FNR > 5) && (NF == 2)){
      split(FILENAME,F,/_/)
      HAPLOTYPES[FILENAME] = HAPLOTYPES[FILENAME] " " substr($2,F[4],1) substr($2,F[5],1)
   }END{
      for (f in HAPLOTYPES) {
         print f "\t" substr(HAPLOTYPES[f],2)
      }
   }' nexus/*.nex > haplotypes
fi

# Some loci include gaps, missing data, recombination events, or more than 2 alleles
# in at least one of the variable sites. I remove them and keep only loci with exactly
# 3 haplotypes, because loci with 2 haplotypes have linked variable sites and loci with
# 4 haplotypes show recombination.

if [ ! -e selected_loci.txt ]; then
   grep -v -P " N|N |-" haplotypes | \
   gawk '{
      for (i=2;i<=NF;i++) {
         FIRST[substr($i,1,1)]=1
         SECOND[substr($i,2,1)]=1
         HAP[$i]=1
      }
      for (f in FIRST) A++
      for (s in SECOND) B++
      for (h in HAP) C++
      delete(FIRST); delete(SECOND); delete(HAP)
      if ((A == 2) && (B==2) && (C == 3)) print $0
      A=0; B=0; C=0
   }' > selected_loci.txt
fi
