BEGIN{
   SEQUENCES = ""
   SPECIES = 0
   INFORMATIVE = 0
   FILENUM = 0
   OFFSET = 0
   NTAX = 0
   NCHAR = 0
}(/^\/\//){
   for (s in S) SPECIES++
   if (SPECIES == 4) {
      INFORMATIVE = gsub(/*/, "*")
      if (INFORMATIVE > 1) {
         FILENUM++
         OUTFILE = sprintf("se%04u_N%u_I%u", FILENUM, SPECIES, INFORMATIVE)
         LINE = $0
         for (i = 1; i <= INFORMATIVE; i++) {
            POSITION = index(LINE, "*") - OFFSET
            OUTFILE = OUTFILE "_" POSITION
            sub(/*/,"x",LINE)
         }
         OUTFILE = OUTFILE ".nex"
         print "#NEXUS" >OUTFILE
         print "BEGIN DATA;" >OUTFILE
         print "DIMENSIONS NTAX=" NTAX " NCHAR=" NCHAR ";" >OUTFILE
         print "FORMAT DATATYPE = DNA GAP = - MISSING = N;" >OUTFILE
         print "MATRIX" >OUTFILE
         print substr(SEQUENCES, 2) >OUTFILE
         print ";" >OUTFILE
         print "END;" >OUTFILE
      }
   }
   SPECIES = 0
   INFORMATIVE = 0
   delete(S)
   NTAX = 0
   NCHAR = 0
   SEQUENCES = ""
}(/^>/){
   S[substr($1,2,4)] = 1
   OFFSET = index($0, " " substr($2, 1, 1))
   NTAX++
   NCHAR = length($2)
   SEQUENCES = SEQUENCES "\n" substr($0,2)
}
