BEGIN{
   SP[1] = "BlCa"
   SP[2] = "BlCa"
   SP[3] = "BlCa"
   SP[4] = "BlCa"
   SP[5] = "BlCa"
   SP[6] = "BlCa"
   SP[7] = "BlCl"
   SP[8] = "BlCl"
   SP[9] = "BlCl"
   SP[10] = "BlCl"
   SP[11] = "BlCl"
   SP[12] = "BlCl"
   SP[13] = "StCa"
   SP[14] = "StCa"
   SP[15] = "StCa"
   SP[16] = "StCa"
   SP[17] = "StCa"
   SP[18] = "StCf"
   SP[19] = "StCf"
   SP[20] = "StCf"
   SP[21] = "StCf"
   SP[22] = "StCf"
   SP[23] = "StCf"
} (((NF == 6) || (NF == 7)) && (/\|/)) {
   LIST[SP[$(NF - 5)] SP[$(NF - 4)] "|" SP[$(NF - 2)] SP[$(NF - 1)]] = LIST[SP[$(NF - 5)] SP[$(NF - 4)] "|" SP[$(NF - 2)] SP[$(NF - 1)]] "," $NF
}END{
   for (Q in LIST) {
      split(LIST[Q],SCORES,",")
      for (i in SCORES) {
         print Q "\t" SCORES[i]
      }
   }
}
