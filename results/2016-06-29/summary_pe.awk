(FNR == 1){
   split(FILENAME,A,/\./)
   SAMPLE = A[1]
}
(/were paired; of these:/){
   TOTAL = $1
}(/aligned concordantly 0 times/){
   NOT_CONCORDANT = $1
}(/aligned concordantly exactly 1 time/){
   CONCORDANT = $1
}(/aligned concordantly >1 times/){
   AMBIGUOUS = $1
}(/aligned discordantly 1 time/){
   DISCORDANT = $1
}(/aligned 0 times/){
   UNMAPPED_MATES = $1
}(/aligned exactly 1 time/){
   MAPPED_MATES = $1
}(/aligned >1 times/){
   AMBIGUOUS_MATES = $1
}(/overall alignment rate/){
   printf "%s\t%9d\t%10d\t%10d\t%8d\t%8d\t%8d\t%10d\t%10d\t%7d\n", \
           SAMPLE, TOTAL, CONCORDANT, AMBIGUOUS, DISCORDANT, NOT_CONCORDANT - DISCORDANT, UNMAPPED_MATES, MAPPED_MATES, AMBIGUOUS_MATES, $1
}
