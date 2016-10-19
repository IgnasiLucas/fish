(FNR == 1){
   split(FILENAME,A,/\./)
   SAMPLE = A[1]
}
(/reads; of these:/){
   TOTAL = $1
}(/aligned 0 times/){
   UNMAPPED = $1
}(/aligned exactly 1 time/){
   UNIQUE = $1
}(/aligned >1 times/){
   AMBIGUOUS = $1
}(/overall alignment rate/){
   printf "%s\t%9d\t%8d\t%9d\t%9d\t%7d\n", SAMPLE, TOTAL, UNMAPPED, UNIQUE, AMBIGUOUS, $1
}
