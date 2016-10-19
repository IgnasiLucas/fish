((/^[^#]/) && (FILENAME ~ /se/)){
   TOTAL_SE += $2
   UNIQ_SE  += $4
   MULT_SE  += $5
   UNMAP_SE += $3
}((/^[^#]/) && (FILENAME ~ /pe/)){
   TOTAL_PE += $2
   CONC_UNI += $3
   CONC_MUL += $4
   DISCORD  += $5
   UNMAP_PE += $6

   REMAIN_MATES += $6 * 2
   UNIQUE_MATES += $8
   MULTI_MATES  += $9
   UNMAP_MATES  += $7
}END{
   printf("% 8u          \n", TOTAL_SE)
   printf("% 8u (% 6.2f \%)\n", UNIQ_SE, 100 * UNIQ_SE/TOTAL_SE)
   printf("% 8u (% 6.2f \%)\n", MULT_SE, 100 * MULT_SE/TOTAL_SE)
   printf("% 8u (% 6.2f \%)\n", UNMAP_SE, 100 * UNMAP_SE/TOTAL_SE)
   print ""
   printf("% 8u          \n", TOTAL_PE)
   printf("% 8u (% 6.2f \%)\n", CONC_UNI, 100 * CONC_UNI/TOTAL_PE)
   printf("% 8u (% 6.2f \%)\n", CONC_MUL, 100 * CONC_MUL/TOTAL_PE)
   printf("% 8u (% 6.2f \%)\n", DISCORD, 100 * DISCORD/TOTAL_PE)
   printf("% 8u (% 6.2f \%)\n", UNMAP_PE, 100 * UNMAP_PE/TOTAL_PE)
   print ""
   printf("% 8u          \n", REMAIN_MATES)
   printf("% 8u (% 6.2f \%)\n", UNIQUE_MATES, 100 * UNIQUE_MATES/REMAIN_MATES)
   printf("% 8u (% 6.2f \%)\n", MULTI_MATES, 100 * MULTI_MATES/REMAIN_MATES)
   printf("% 8u (% 6.2f \%)\n", UNMAP_MATES, 100 * UNMAP_MATES/REMAIN_MATES)
}
