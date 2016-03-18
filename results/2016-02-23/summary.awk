BEGIN{
   i = 0
}(FNR == 1) {
   i++
   ClustTh = substr(FILENAME, 21, 2)
   THRESHOLD[i] = ClustTh
}((FNR <= 26) && (/^StC|^BlC/)){
   SAMPLES[$1] = 1
   ClustNum[ClustTh "," $1] = $2
   MeanDepth[ClustTh "," $1] = $3
   ClustNumD3[ClustTh "," $1] = $5
   MeanDepthD3[ClustTh "," $1] = $6
}END{
   HEADER = ""
   for (j = 1; j <= i; j++){
      HEADER = HEADER "\tClustNum_" THRESHOLD[j] "\tMeanDepth_" THRESHOLD[j] "\tClustNumD3_" THRESHOLD[j] "\tMeanDepthD3_" THRESHOLD[j]
   }
   print HEADER
   for (s in SAMPLES) {
      LINE = s
      for (j = 1; j <= i; j++){
         LINE = LINE "\t" ClustNum[THRESHOLD[j] "," s] "\t" MeanDepth[THRESHOLD[j] "," s] "\t" ClustNumD3[THRESHOLD[j] "," s] "\t" MeanDepthD3[THRESHOLD[j] "," s]
      }
      print LINE
   }
}
