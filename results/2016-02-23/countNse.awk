# Counts proportions of Ns in each merged read.
(/^>/){
   if (L > 0) {
      printf("%.8f\n", N/L)
   }
   L = 0
   N = 0
}(/^[^>]/){
   L += length($1)
   N += gsub(/N/,"n")
}END{
   printf("%.8f\n", N/L)
}
