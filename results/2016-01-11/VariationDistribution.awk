BEGIN{
   N = 0
   LENGTH = 0
   for (i = 1; i <= 24; i++) {
      TOTALPOSITIONS["read1", i] = 0
      TOTALPOSITIONS["overlap", i] = 0
      TOTALPOSITIONS["read2", i] = 0
      VARPOSITIONS["read1", i] = 0
      VARPOSITIONS["overlap", i] = 0
      VARPOSITIONS["read2", i] = 0
   }
}(/^\/\//){
   LINE = $0
   VARIANT = gsub(/*|-/,"x", LINE)
   if (VARIANT > 0) {
      if (LENGTH <= 300) {
         TOTALPOSITIONS["overlap", N] += LENGTH
         VARPOSITIONS["overlap", N] += VARIANT
      } else {
         for (i = 1; i <= VARIANT; i++) {
            VARPOS = index(LINE, "x") - OFFSET
            if (VARPOS <= LENGTH - 300) {
               VARPOSITIONS["read1", N]++
            } else {
               if (VARPOS <= 300) {
                  VARPOSITIONS["overlap", N]++
               } else {
                  VARPOSITIONS["read2", N]++
               }
            }
            sub(/x/,"-",LINE)
         }
         TOTALPOSITIONS["read1", N] += LENGTH - 300
         TOTALPOSITIONS["overlap", N] += 600 - LENGTH
         TOTALPOSITIONS["read2", N] += LENGTH - 300
      }
   } else {
      if (LENGTH <= 300) {
         TOTALPOSITIONS["overlap", N] += LENGTH
      } else {
         TOTALPOSITIONS["read1", N] += LENGTH - 300
         TOTALPOSITIONS["overlap", N] += 600 - LENGTH
         TOTALPOSITIONS["read2", N] += LENGTH - 300
      }
   }
   N = 0
   LENGTH = 0
   delete(POSITIONS)
}(/^>/){
   N++
   LENGTH = length($2)
   OFFSET = index($0, " " substr($2,1,1))
}END{
   print "\tRead1\t\tOverlap\t\tRead2"
   print "Sequences\tVariable\tTotal\tVariable\tTotal\tVariable\tTotal"
   for (i = 1; i <= 24; i++) {
      printf("%u\t%u\t%u\t%u\t%u\t%u\t%u\n", \
             i,
             VARPOSITIONS["read1", i], TOTALPOSITIONS["read1", i],
             VARPOSITIONS["overlap", i], TOTALPOSITIONS["overlap", i],
             VARPOSITIONS["read2", i], TOTALPOSITIONS["read2", i])
   }
}
