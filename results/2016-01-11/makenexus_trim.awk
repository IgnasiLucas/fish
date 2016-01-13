# This versions will use only sequences MINLENGTH bp or longer, and it will
# substitute bases between the STARTRM and ENDRM (floats) of their length by
# Ns. It will remove the variable sites thought to be there from the name of
# the file.
BEGIN{
   SEQUENCES = ""
   SPECIES = 0
   INFORMATIVE = 0
   FILENUM = 0
   OFFSET = 0
   NTAX = 0
   NCHAR = 0
   if (MINLENGTH == "") MINLENGTH = 0
   if (STARTRM == "") STARTRM = 1.1
   if (ENDRM == "") ENDRM = 1.1
}(/^\/\//){
   for (s in S) SPECIES++
   if ((SPECIES == 4) && (NCHAR >= MINLENGTH)) {
      INFORMATIVE = gsub(/*/, "*")
      RM1 = sprintf("%.0f", NCHAR * STARTRM)
      RM2 = sprintf("%.0f", NCHAR * ENDRM)
      if (INFORMATIVE > 1) {
         LINE = $0
         NEWINFORMATIVE = INFORMATIVE
         for (i = 1; i <= INFORMATIVE; i++) {
            POSITION = index(LINE, "*") - OFFSET
            if ((POSITION > RM1) && (POSITION <= RM2)) {
               NEWINFORMATIVE--
            } else {
               SUFIX = SUFIX "_" POSITION
            }
            sub(/*/,"x",LINE)
         }
      }
      if (NEWINFORMATIVE > 1) {
         FILENUM++
         OUTFILE = sprintf("se%04u_N%u_I%u%s.nex", FILENUM, SPECIES, NEWINFORMATIVE, SUFIX)
         print "#NEXUS" >OUTFILE
         print "BEGIN DATA;" >OUTFILE
         print "DIMENSIONS NTAX=" NTAX " NCHAR=" NCHAR ";" >OUTFILE
         print "FORMAT DATATYPE = DNA GAP = - MISSING = N;" >OUTFILE
         print "MATRIX" >OUTFILE
         split(substr(SEQUENCES, 2), ALIGNMENT, "\n")
         for (NAMESEQ in ALIGNMENT) {
            split(ALIGNMENT[NAMESEQ], SEQ)
            NEWSEQ = substr(SEQ[2], 1, RM1)
            for (i = RM1 + 1; i <= RM2; i++) NEWSEQ = NEWSEQ "N"
            NEWSEQ = NEWSEQ substr(SEQ[2], RM2 + 1)
            print SEQ[1] "  " NEWSEQ >OUTFILE
         }
         print ";" >OUTFILE
         print "END;" >OUTFILE
      }
   }
   SPECIES = 0
   INFORMATIVE = 0
   NEWINFORMATIVE = 0
   delete(S)
   NTAX = 0
   NCHAR = 0
   SEQUENCES = ""
   SUFIX = ""
}(/^>/){
   S[substr($1,2,4)] = 1
   OFFSET = index($0, " " substr($2, 1, 1))
   NTAX++
   NCHAR = length($2)
   SEQUENCES = SEQUENCES "\n" substr($0,2)
}
