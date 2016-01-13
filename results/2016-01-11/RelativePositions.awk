BEGIN{
   if (MAX == 0) MAX = 100000
}(/^\/\//){
   if ((LENGTH >= MIN) && (LENGTH <= MAX)) {
      LINE = $0
      VARIANT = gsub(/*|-/, "x", LINE)
      for (i = 1; i <= VARIANT; i++) {
         print (index(LINE, "x") - OFFSET) / LENGTH
         sub(/x/, "-", LINE)
      }
   }
}(/^>/){
   OFFSET = index($0, " " substr($2, 1, 1))
   LENGTH = length($2)
}
