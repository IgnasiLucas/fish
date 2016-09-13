#!/bin/bash
#
#				2016-06-29
#				----------
#
# The consensus sequences of the read pairs that could be merged are ready in
# 2016-05-10/clust.86/pooled.consens.gz. It's a fasta file. Now, the idea is to
# use it as a reference genome, and map to it all the reads. Ideally, every
# consensus sequence should be covered by at least one merged read. The success
# of the approach can be measured by the proportion of non-merged reads that also
# map (unambiguously) to the consensus sequences, because these would be reads
# that failed to merge due to errors and could not be clustered with their
# homologous reads before.
#
# The fastq files, merged and trimmed, are in 2016-05-03.

DATADIR=`pwd | sed 's/results/data/'`
FASTQDIR=`pwd | sed 's/2016-06-29/2016-05-03/'`
if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi

if [ ! -e $DATADIR/consensus.fa ]; then
   cp ../2016-05-10/clust.86/pooled.consens.gz $DATADIR/
   gunzip $DATADIR/pooled.consens.gz
   mv     $DATADIR/pooled.consens $DATADIR/consensus.fa
fi

SAMPLE=(ZeroNotUsed St0001 St0003        St0015 St0016 St0019 St0037 St0039
                    St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
                    Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

for i in `seq 1 23`; do
   if [ ! -e $DATADIR/${SAMPLE[$i]}.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed.fastq $DATADIR/${SAMPLE[$i]}.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R1.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R1.fastq $DATADIR/${SAMPLE[$i]}_R1.fastq
   fi
   if [ ! -e $DATADIR/${SAMPLE[$i]}_R2.fastq ]; then
      ln -s $FASTQDIR/${SAMPLE[$i]}_trimmed_R2.fastq $DATADIR/${SAMPLE[$i]}_R2.fastq
   fi
done

if [ ! -e fish.1.bt2 ]; then
   bowtie2-build $DATADIR/consensus.fa fish
fi

function map_single_ends {
   bowtie2 --very-sensitive \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1 \
           --un $1.se.unmapped.fastq \
           -x fish -U $DATADIR/$1.fastq \
           -S $1.se.sam &> $1.se.bowtie2.log
}

function map_paired_ends {
   bowtie2 --very-sensitive \
           --maxins 600 \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1} \
           --un-conc $1.pe.unmapped_R%.fastq \
           -x fish -1 $DATADIR/$1'_R1.fastq' -2 $DATADIR/$1'_R2.fastq' \
           -S $1.pe.sam &> $1.pe.bowtie2.log
}

PROC=`grep -P '^processor' /proc/cpuinfo | wc -l`

if [ $PROC -gt 23 ]; then
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]} &
      fi
   done
   wait
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]} &
      fi
   done
   wait
elif [ $PROC -ge 6 ]; then
   for j in 1 7 13; do
      for i in `seq $j $(( j + 5 ))`; do
         if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
            map_single_ends ${SAMPLE[$i]} &
         fi
      done
      wait
      for i in `seq $j $(( j + 5 ))`; do
         if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
            map_paired_ends ${SAMPLE[$i]} &
         fi
      done
      wait
   done
   for i in 19 20 21 22 23; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]} &
      fi
   done
   wait
   for i in 19 20 21 22 23; do
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]} &
      fi
   done
   wait
else
   for i in `seq 1 23`; do
      if [ ! -e ${SAMPLE[$i]}.se.sam ]; then
         map_single_ends ${SAMPLE[$i]}
      fi
      if [ ! -e ${SAMPLE[$i]}.pe.sam ]; then
         map_paired_ends ${SAMPLE[$i]}
      fi
   done
fi

#if [ ! -e summary_se.txt ]; then
   echo "## Tot.Reads: Total number of reads."                             > summary_se.txt
   echo "## Unmapped:  Aligned 0 times."                                  >> summary_se.txt
   echo "## Uniq.Map.: Aligned exactly 1 time."                           >> summary_se.txt
   echo "## Mult.Map.: Aligned >1 times."                                 >> summary_se.txt
   echo "## Overall: Overall mapping rate."                               >> summary_se.txt
   echo -e "# Sample\tTot.Reads\tUnmapped\tUniq.Map.\tMult.Map.\tOverall" >> summary_se.txt
   gawk -f summary_se.awk *.se.bowtie2.log | sort -nrk 2                  >> summary_se.txt
#fi
#if [ ! -e summary_pe.txt ]; then
   echo "## Tot.Pairs: Total number of read pairs."                        > summary_pe.txt
   echo "## Conc.Uniq.: Uniquely and concordantly mapped pairs."          >> summary_pe.txt
   echo "## Conc.Mult.: Multiply and concordantly mapped pairs."          >> summary_pe.txt
   echo "## Discord.:  Discordantly mapped pairs."                        >> summary_pe.txt
   echo "## Un.Pairs: Unmapped pairs."                                    >> summary_pe.txt
   echo "## Un.Mates: Reads or mates from unmapped pairs that do not map either as single reads." >> summary_pe.txt
   echo "## Uniq.Mates: Mates from unmapped pairs that map uniquely as single reads."             >> summary_pe.txt
   echo "## Mult.Mates: Mates from unmapped pairs that map ambiguously as single reads."          >> summary_pe.txt
   echo "## Overall: Overall mapping rate."                               >> summary_pe.txt
   echo -e "# Sample\tTot.Pairs\tConc.Uniq.\tConc.Mult.\tDiscord.\tUn.Pairs\tUn.Mates\tUniq.Mates\tMult.Mates\tOverall" >> summary_pe.txt
   gawk -f summary_pe.awk *.pe.bowtie2.log | sort -nrk 2                  >> summary_pe.txt
#fi
gnuplot < summary_se.gnuplot
gnuplot < summary_pe.gnuplot

# CONCLUSIONS
# -----------
#
# The plots show that most reads that were assembled (pe) can be mapped back
# to the contigs that were made with them. There is some ambiguity, though,
# since ~50% of them map in multiple contigs. Without knowing how well they
# map to the alternative loci, it is difficult to make a statement about
# the adequacy of the similarity threshold used in the clustering step (86%).
# In principle, it is a quite permissive threshold, which should not split the
# same locus in different clusters. Instead, I think we are seeing that bowtie2
# is using very sensitive settings, and finding homology among loci. This is
# likely, since the genome of an ancestor of Salmonids suffered a whole
# duplication some 90 mya (Berthelot et al. 2014).
#
# The problem is with the non-assembled reads, that is, the ones that are kept
# as paired ends. I expected many of them were not assembled just because of
# the lower quality of their 3' ends. Figure summary_pe.png shows that indeed
# there is a number of paired ends that map well and uniquely to the contigs
# made from the merged reads. However, they are less than 20%. Hundreds of
# thousands of reads in each sample still need to be analysed separately with
# pyrad.
#
# Next steps: with the reads that are well mapped (with a minimum quality, and
# whether merged or not), I will use freebayes to get the vcf and analyse the
# variation; with the paired ends that did not map, I will try the new ipyrad.
# It will be challenging to merge the information from both sets.
#
# I refuse to force the paired ends to map better to the artificial contigs
# (e.g. by trimming them). However, it may be worth checking how much better
# a merged read maps to the first than to its second hit, in order to determine
# a convenient threshold of mapping quality. I don't want to be stricter than
# I need. Here, bowtie2 reported only the best mapping of each read. I can take
# a sample of reads with low quality mappings and map them again while asking
# for all the mappings to be shown. For some reason, bowtie2 warns that mapping
# qualities will not be informative in that case, though.
