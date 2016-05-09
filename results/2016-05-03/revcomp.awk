{
   REV = ""
   for (i = length($1); i > 0; i--) {
      REV = REV substr($1,i,1)
   }
   gsub(/A/,"t",REV)
   gsub(/C/,"g",REV)
   gsub(/G/,"c",REV)
   gsub(/T/,"a",REV)
   print toupper(REV)
}
