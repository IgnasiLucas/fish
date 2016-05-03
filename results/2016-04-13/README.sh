#!/bin/bash
#
#				2016-04-13
#				----------
#
# It is clear that some, if not most, non-merged read pairs failed to
# merge not because they come from long templates, but because of the
# poor quality of the bases near the end of the read.
#
# Before trying to combine merged and non-merged reads in the same run
# of pyrad, which would be messy, I attempted to merge again with
# more sensitive settings. I used bbmerge, from the bbmap package, but
# it is not as sensitive as pear (although it may be faster). Trying
# again with PEAR, with one sample, I could not increase much the merging
# rate, and convinced myself that I need another solution
#
# The bbmap suite includes tadpole, which build assemblies, quite fast,
# and in a very conservative way. One option is to use the raw reads from
# maybe all the samples to build an assembly. I expect contigs to be either
# single reads that do not merge, or merged pairs of reads. Then, individual
# reads could be mapped to the contigs.
#

FASTQ=$(dirname $(pwd))"/2015-12-16b/joined"
SAMPLE=(ZeroNotUsed St0001 St0003 St0006 St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)
BARCODE=(ZeroNotUsed AGCTA TCGAG TAACCTG ACCGAGT TGGGTGCC GTCTTGCG
         CAATCC GTTCAA AGCTA TCGAG TAACCTG ACCGAGT
         TGGGTGCC GTCTTGCG CAATCC GTTCAA AGCTA TCGAG
         TAACCTG ACCGAGT TGGGTGCC GTCTTGCG CAATCC GTTCAA)
i=19

# I switch to a checkpoints-based conditional execution, to make it easier
# to erase files after being used.

touch checkpoints
for j in 1 2; do
   if ! grep -q ${SAMPLE[$i]}_R$j.fastq checkpoints; then
      echo ${SAMPLE[$i]}_R$j.fastq >> checkpoints
      cp $FASTQ/${SAMPLE[$i]}_R$j.fastq.gz ${SAMPLE[$i]}_R$j.fastq.gz
      gunzip ${SAMPLE[$i]}_R$j.fastq.gz
   fi
done

if ! grep -q ${SAMPLE[$i]}_demultiplexed_R1.fastq checkpoints; then
   echo ${SAMPLE[$i]}_demultiplexed_R1.fastq >> checkpoints
   echo "${BARCODE[$i]} ${SAMPLE[$i]}_demultiplexed_R1.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq" > ${SAMPLE[$i]}.barcode
   sabre pe -m 1 \
            -f ${SAMPLE[$i]}_R1.fastq \
            -r ${SAMPLE[$i]}_R2.fastq \
            -b ${SAMPLE[$i]}.barcode \
            -u ${SAMPLE[$i]}_unknown_R1.fastq \
            -w ${SAMPLE[$i]}_unknown_R2.fastq > ${SAMPLE[$i]}_sabre.log
   echo ${BARCODE[$i]} | gawk -v S=${SAMPLE[$i]} '{print substr($1,2) " " S "_emultiplexed_R1.fastq " S "_emultiplexed_R2.fastq"}' > ${SAMPLE[$i]}.arcode
   sabre pe -m 0 \
            -f ${SAMPLE[$i]}_unknown_R1.fastq \
            -r ${SAMPLE[$i]}_unknown_R2.fastq \
            -b ${SAMPLE[$i]}.arcode \
            -u ${SAMPLE[$i]}_uunknown_R1.fastq \
            -w ${SAMPLE[$i]}_uunknown_R2.fastq > ${SAMPLE[$i]}_ssabre.log
   if [ -s ${SAMPLE[$i]}_emultiplexed_R1.fastq ]; then
      cat_seqs -o ${SAMPLE[$i]}_dem_R1.fastq ${SAMPLE[$i]}_demultiplexed_R1.fastq ${SAMPLE[$i]}_emultiplexed_R1.fastq
      cat_seqs -o ${SAMPLE[$i]}_dem_R2.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq ${SAMPLE[$i]}_emultiplexed_R2.fastq
      mv ${SAMPLE[$i]}_dem_R1.fastq ${SAMPLE[$i]}_demultiplexed_R1.fastq
      mv ${SAMPLE[$i]}_dem_R2.fastq ${SAMPLE[$i]}_demultiplexed_R2.fastq
   fi
   rm ${SAMPLE[$i]}_emultiplexed_*
   rm ${SAMPLE[$i]}_R1.fastq
   rm ${SAMPLE[$i]}_R2.fastq
   rm ${SAMPLE[$i]}_unknown*
   rm ${SAMPLE[$i]}_uunknown*
fi

if ! grep -q ${SAMPLE[$i]}_trimmed_R1.fastq checkpoints; then
   echo ${SAMPLE[$i]}_trimmed_R1.fastq >> checkpoints
   BARCODEREV=`echo ${BARCODE[$i]} | gawk -f revcomp.awk`
   cutadapt -a CGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A "TA"$BARCODEREV"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
            -o ${SAMPLE[$i]}_trimmed_R1.fastq \
            -p ${SAMPLE[$i]}_trimmed_R2.fastq \
            ${SAMPLE[$i]}_demultiplexed_R1.fastq \
            ${SAMPLE[$i]}_demultiplexed_R2.fastq > cutadapt.log
   if [ -s ${SAMPLE[$i]}_trimmed_R1.fastq ]; then
      rm ${SAMPLE[$i]}_demultiplexed*
   fi
fi

if ! grep -q ${SAMPLE[$i]}_interleaved.fastq checkpoints; then
   reformat.sh in=${SAMPLE[$i]}_trimmed_R1.fastq \
               in2=${SAMPLE[$i]}_trimmed_R2.fastq \
               out=${SAMPLE[$i]}_interleaved.fastq \
               verifypaired=t \
               minlength=100 \
               rcompmate=t
fi

if [ ! -e ${SAMPLE[$i]}.fa ]; then
   tadpole.sh in=${SAMPLE[$i]}_interleaved.fastq \
              out=${SAMPLE[$i]}.fa \
              k=25 \
              mincountseed=1 \
              mincountextend=1 \
              branchmult1=1 \
              branchmult2=1
fi

# Conclusion
# ----------
#
# k-mer based assembly tools will not work well with this type of data, where some read pairs
# can be merged and other cannot, with errors preventing some pairs from merging. It is better
# to merge the pairs that can be merged with PEAR, and if anything, build a reference with them
# and map the rest.
