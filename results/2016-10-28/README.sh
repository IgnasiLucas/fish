#!/bin/bash
#
#				2016-10-28
#				----------
#
# Here I abandon the idea that a reference set of loci can be built from
# the pooled merged reads, and I use the previous efforts to distinguish the
# paired end (non-merged) reads that can be safely used in parallel to the
# merged reads, without risk of duplicating loci. That is, only paired end
# reads that do not show any detectable homology with merged reads should
# be assumed to come from different loci. I will therefore run ipyrad twice
# independently, once for merged reads and once for non-merged-and-non-mappable
# reads.
#

  DATADIR=`pwd | sed 's/results/data/'`
PAIREDDIR=`pwd | sed 's/2016-10-28/2016-10-19/'`
MERGEDDIR=`pwd | sed 's|2016-10-28|2016-10-03/merged|'`

if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi
if [ ! -d $DATADIR/merged ]; then mkdir $DATADIR/merged; fi
if [ ! -d $DATADIR/paired ]; then mkdir $DATADIR/paired; fi

SAMPLE=(ZeroNotUsed StCa0001 StCa0003          StCa0015 StCa0016 StCa0019 StCf0037 StCf0039
                    StCf0043 StCf0044 StCf0049 StCf0050 BlCa0065 BlCa0076 BlCa0080 BlCa0083
                    BlCa0104 BlCa0108 BlCl0091 BlCl0093 BlCl0094 BlCl0095 BlCl0098 BlCl0116)

for i in `seq 1 23`; do
   if [ ! -e $DATADIR/merged/${SAMPLE[$i]}.fastq ]; then
      ln -s $MERGEDDIR/${SAMPLE[$i]}.assembled.fastq $DATADIR/merged/${SAMPLE[$i]}.fastq
   fi
   if [ ! -e $DATADIR/paired/${SAMPLE[$i]}_R1_.fastq ]; then
      ln -s $PAIREDDIR/${SAMPLE[$i]}.pe.unmapped_R1.fastq $DATADIR/paired/${SAMPLE[$i]}_R1_.fastq
   fi
   if [ ! -e $DATADIR/paired/${SAMPLE[$i]}_R2_.fastq ]; then
      ln -s $PAIREDDIR/${SAMPLE[$i]}.pe.unmapped_R2.fastq $DATADIR/paired/${SAMPLE[$i]}_R2_.fastq
   fi
done

if [ ! -e populations.txt ]; then
   for i in `seq 1 23`; do
      echo "${SAMPLE[$i]} ${SAMPLE[$i]:0:4}" >> populations.txt
   done
fi

if [ ! -e params-merged.txt ]; then
   ipyrad -n merged
   # This only adds the end-of-line before the end of file
   echo >> params-merged.txt
#                    ------- ipyrad params file (v.0.4.1)--------------------------------------------
   sed -i "/## \[1\]/c merged                         ## [1] [project_dir]: Project dir (made in curdir if not present)"           params-merged.txt
   sed -i "/## \[4\]/c $DATADIR/merged/*.fastq        ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files"    params-merged.txt
   sed -i "/## \[7\]/c ddrad                          ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc."               params-merged.txt
   sed -i "/## \[8\]/c TAA,CGG                        ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)" params-merged.txt
   sed -i "/## \[9\]/c 5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read"    params-merged.txt
   sed -i "/## \[10\]/c 30                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and...)" params-merged.txt
   sed -i "/## \[11\]/c 6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling"     params-merged.txt
   sed -i "/## \[12\]/c 1                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling"       params-merged.txt
   sed -i "/## \[13\]/c 100                            ## [13] [maxdepth]: Max cluster depth within samples"                       params-merged.txt
   sed -i "/## \[14\]/c 0.86                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly"       params-merged.txt
   sed -i "/## \[16\]/c 0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)"   params-merged.txt
   sed -i "/## \[17\]/c 35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim"      params-merged.txt
   sed -i "/## \[18\]/c 2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences" params-merged.txt
   sed -i "/## \[19\]/c 5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)"   params-merged.txt
   sed -i "/## \[20\]/c 88, 88                         ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)"     params-merged.txt
   sed -i "/## \[21\]/c 1                              ## [21] [min_samples_locus]: Min # samples per locus for output"            params-merged.txt
   sed -i "/## \[22\]/c 40, 40                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)"                    params-merged.txt
   sed -i "/## \[23\]/c 8, 8                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)"             params-merged.txt
   sed -i "/## \[24\]/c 0.9                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)" params-merged.txt
   sed -i "/## \[25\]/c 0, 0                           ## [25] [edit_cutsites]: Edit cut-sites (R1, R2) (see docs)"                params-merged.txt
   sed -i "/## \[26\]/c 4, 4, 4, 4                     ## [26] [trim_overhang]: Trim overhang (see docs) (R1>, <R1, R2>, <R2)"     params-merged.txt
   sed -i "/## \[27\]/c *                              ## [27] [output_formats]: Output formats (see docs)"                        params-merged.txt
   sed -i "/## \[28\]/c populations.txt                ## [28] [pop_assign_file]: Path to population assignment file"              params-merged.txt
fi

if [ ! -e params-paired.txt ]; then
   ipyrad -n paired
   echo >> params-paired.txt
#                    ------- ipyrad params file (v.0.4.1)--------------------------------------------
   sed -i "/## \[1\]/c paired                         ## [1] [project_dir]: Project dir (made in curdir if not present)"           params-paired.txt
   sed -i "/## \[4\]/c $DATADIR/paired/*.fastq        ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files"    params-paired.txt
   sed -i "/## \[7\]/c pairddrad                      ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc."               params-paired.txt
   sed -i "/## \[8\]/c TAA,CGG                        ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)" params-paired.txt
   sed -i "/## \[9\]/c 5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read"    params-paired.txt
   sed -i "/## \[10\]/c 30                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and...)" params-paired.txt
   sed -i "/## \[11\]/c 6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling"     params-paired.txt
   sed -i "/## \[12\]/c 1                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling"       params-paired.txt
   sed -i "/## \[13\]/c 100                            ## [13] [maxdepth]: Max cluster depth within samples"                       params-paired.txt
   sed -i "/## \[14\]/c 0.86                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly"       params-paired.txt
   sed -i "/## \[16\]/c 0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)"   params-paired.txt
   sed -i "/## \[17\]/c 35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim"      params-paired.txt
   sed -i "/## \[18\]/c 2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences" params-paired.txt
   sed -i "/## \[19\]/c 5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)"   params-paired.txt
   sed -i "/## \[20\]/c 88, 88                         ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)"     params-paired.txt
   sed -i "/## \[21\]/c 1                              ## [21] [min_samples_locus]: Min # samples per locus for output"            params-paired.txt
   sed -i "/## \[22\]/c 40, 40                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)"                    params-paired.txt
   sed -i "/## \[23\]/c 8, 8                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)"             params-paired.txt
   sed -i "/## \[24\]/c 0.9                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)" params-paired.txt
   sed -i "/## \[25\]/c 0, 0                           ## [25] [edit_cutsites]: Edit cut-sites (R1, R2) (see docs)"                params-paired.txt
   sed -i "/## \[26\]/c 4, 4, 4, 4                     ## [26] [trim_overhang]: Trim overhang (see docs) (R1>, <R1, R2>, <R2)"     params-paired.txt
   sed -i "/## \[27\]/c *                              ## [27] [output_formats]: Output formats (see docs)"                        params-paired.txt
   sed -i "/## \[28\]/c populations.txt                ## [28] [pop_assign_file]: Path to population assignment file"              params-paired.txt
fi

for STEP in 2 3 4 5 6 7; do
   if ipyrad -p params-merged.txt -r | gawk -v STEP=$STEP '((/^step/) && ($2 == STEP ":")){print}' | tee z1 | grep -q None; then
      ipyrad -p params-merged.txt -s $STEP -c 24 1> merged/step$STEP.log 2> merged/step$STEP.err
      # This is to reduce the size of the log:
      gawk -v FS="\r" '{print $NF}' merged/step$STEP.log > z2
      mv z2 merged/step$STEP.log
   else
      cat z1
      rm z1
   fi

   if ipyrad -p params-paired.txt -r | gawk -v STEP=$STEP '((/^step/) && ($2 == STEP ":")){print}' | tee z1 | grep -q None; then
      ipyrad -p params-paired.txt -s $STEP -c 24 1> paired/step$STEP.log 2> paired/step$STEP.err
      gawk -v FS="\r" '{print $NF}' paired/step$STEP.log > z2
      mv z2 paired/step$STEP.log
   else
      cat z1
      rm z1
   fi
done
rm z1
