set terminal png large
set output 'summary_pe.png'
set style data histograms
set style histogram rowstacked
set style fill solid 1.0 border -1
set ylabel "Number of reads"
set bmargin 5
set xtics nomirror rotate by 45 right

plot 'summary_pe.txt' using 6 t "Unmapped", '' using 3 t "Concordant", '' using 4 t "Multiple", '' using 5:xtic(1) t "Discordant"

