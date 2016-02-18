#!/bin/bash
#				2016-02-18
#				----------
#
# I will use popart (Leigh & Bryant 2015) to represent at least some haplotype
# networks. It requires the nexus files to be in different format, with only the
# unique haplotypes in the data block and then the number of times each haplotype
# appears in each population.

NEXUSDIR=`pwd | sed 's_2016-02-18_2016-01-13/nexus_'`

if [ ! -d nexus ]; then
   mkdir nexus
fi

for i in `ls -1 $NEXUSDIR`; do
   if [ ! -e nexus/$i ]; then
      gawk '((NR > 5) && ($0 !~ /;/)){
         SPECIES = substr($1,1,4)
         NCHAR = length($2)
         FREQUENCY[SPECIES "," $2]++
         HAPLOTYPES[$2]++
         SPLIST[SPECIES]++
      }END{
         for (h in HAPLOTYPES) NTAX++
         print "#NEXUS"
         print "BEGIN DATA;"
         print "   DIMENSIONS NTAX=" NTAX " NCHAR=" NCHAR ";"
         print "   FORMAT DATATYPE=DNA MISSING=N GAP=-;"
         print "   MATRIX"
         SEQID=1
         for (h in HAPLOTYPES) {
            printf("seq_%02u      %s\n", SEQID, h)
            ID2SEQ[SEQID] = h
            SEQID++
         }
         print ";"
         print ""
         print "BEGIN TRAITS;"
         print "   Dimensions NTRAITS=4;"
         print "   Format labels=yes missing=? separator=Comma;"
         print "   TraitLabels BlCa BlCl StCa StCf;"
         print "   Matrix"
         for (i = 1; i < SEQID; i++) {
            printf("seq_%02u %u,%u,%u,%u\n", i, FREQUENCY["BlCa," ID2SEQ[i]] + 0, \
                                                FREQUENCY["BlCl," ID2SEQ[i]] + 0, \
                                                FREQUENCY["StCa," ID2SEQ[i]] + 0, \
                                                FREQUENCY["StCf," ID2SEQ[i]] + 0)
         }
         print ";"
         print "END;"
      }' $NEXUSDIR/$i > nexus/$i
   fi
done
