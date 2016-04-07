#!/bin/bash
#
# Sibelle Vilaça asked me if all the reads in a .derep file are present in the
# corresponding .clustS file. She noticed some reads missing, and she's using
# merged reads. Let's count reads in those pairs of files, in the se folder.


SAMPLE=(ZeroNotUsed StCa0001 StCa0003 StCa0006 StCa0015 StCa0016 StCa0019 StCf0037 StCf0039
                    StCf0043 StCf0044 StCf0049 StCf0050 BlCa0065 BlCa0076 BlCa0080 BlCa0083
                    BlCa0104 BlCa0108 BlCl0091 BlCl0093 BlCl0094 BlCl0095 BlCl0098 BlCl0116)

if [ ! -e CheckMissingReads.txt ]; then
   echo -e "#Sample\tClustering\tReadsDerep\tReadsClustS\tDifference" > CheckMissingReads.txt
   for i in `seq 1 24`; do
      DEREP=`grep ">" se/edits/${SAMPLE[$i]}.derep | wc -l`
      for j in 74 76 78 80 82 84 86 88 90 92 94 96 98; do
         CLUSTS=`zless se/clust.$j/${SAMPLE[$i]}.clustS.gz | grep ">" | wc -l`
         echo -e "${SAMPLE[$i]}\t.$j\t$DEREP\t$CLUSTS\t$(( $DEREP - $CLUSTS ))" >> CheckMissingReads.txt
      done
   done
fi
