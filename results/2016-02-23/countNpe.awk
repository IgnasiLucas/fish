# Prints out two columns, with proportions of Ns in first read
# and proportion of Ns in both reads.
(/^>/){
   if (S != "") {
      split(S, A, /nnnn/)
      L1 = length(A[1])
      L2 = L1 + length(A[2])
      N1 = gsub(/N/, "n", A[1])
      N2 = N1 + gsub(/N/, "n", A[2])
      for (i = 1; i <= R; i++){
         printf("%.8f\t%.8f\n", N1/L1, N2/L2)
      }
   }
   S = ""
   split($1, A, /;|=/)
   R = A[3]
}(/^[^>]/){
   S = S $1
}END{
   split(S, A, /nnnn/)
   L1 = length(A[1])
   L2 = L1 + length(A[2])
   N1 = gsub(/N/, "n", A[1])
   N2 = N1 + gsub(/N/, "n", A[2])
   for (i = 1; i <= R; i++){
      printf("%.8f\t%.8f", N1/L1, N2/L2)
   }
}
