#!/bin/bash
#
#				2015-06-15
#				----------
#
# I am digesting in silico the genome sequence of Salmo salar as a proxy
# to the Coregonus genome.

DATADIR=`pwd | sed "s/results/data/"`
if [ ! -d $DATADIR ]; then
	mkdir $DATADIR
fi

if [ ! -e $DATADIR/Salmo_salar.fna ]; then
	wget ftp://ftp.ncbi.nih.gov/genomes/genbank/vertebrate_other/Salmo_salar/latest_assembly_versions/GCA_000233375.3_ICSASG_v1/GCA_000233375.3_ICSASG_v1_genomic.fna.gz
	gunzip GCA_000233375.3_ICSASG_v1_genomic.fna.gz
	mv GCA_000233375.3_ICSASG_v1_genomic.fna $DATADIR/Salmo_salar.fna
fi

touch checkpoints

# I tried to use the R package SimRAD, but it requires too much memory and
# does not work.
if ! grep -q SimRAD checkpoints; then
	echo SimRAD >> checkpoints
	R --save < plot.R
fi

# I'll try now with Marie's python program.
if ! grep -q wide checkpoints; then
	echo wide >> checkpoints
	echo -e "MseI\tT^TA_A\tMspI\tC^CG_G" > enzymes
	python ~/Software/digestion/digestion/RAD_digestion.py \
		-e enzymes \
		-g $DATADIR/Salmo_salar.fna \
		-o wide/ \
		--dd \
		--min 450 \
		--max 750
fi

if ! grep -q narrow_500-550 checkpoints; then
	echo narrow_500-550 >> checkpoints
	python ~/Software/digestion/digestion/RAD_digestion.py \
		-e enzymes \
		-g $DATADIR/Salmo_salar.fna \
		-o narrow_500-550/ \
		--dd \
		--min 500 \
		--max 550
fi

if ! grep -q narrow_550-600 checkpoints; then
        echo narrow_550-600 >> checkpoints
        python ~/Software/digestion/digestion/RAD_digestion.py \
                -e enzymes \
                -g $DATADIR/Salmo_salar.fna \
                -o narrow_550-600/ \
                --dd \
                --min 550 \
                --max 600
fi

if ! grep -q narrow_600-650 checkpoints; then
        echo narrow_600-650 >> checkpoints
        python ~/Software/digestion/digestion/RAD_digestion.py \
                -e enzymes \
                -g $DATADIR/Salmo_salar.fna \
                -o narrow_600-650/ \
                --dd \
                --min 600 \
                --max 650
fi

if ! grep -q all_sizes checkpoints; then
	echo all_sizes >> checkpoints
	python ~/Software/digestion/digestion/RAD_digestion.py \
		-e enzymes \
		-g $DATADIR/Salmo_salar.fna \
		-o all_sizes/ \
		--min 1 \
		--max 1000000
fi

if [ ! -e wide/lengths_inner.hist ]; then
	grep ">" wide/MseI_MspI.fasta | \
	gawk '($1 !~ /End/){split($1,A,/[-\|]/); F[A[4]-A[3]]++}END{for (f in F) print f "\t" F[f]}' | \
	sort -nk 1 > wide/lengths_inner.hist
fi

if [ ! -e wide/lengths_all.hist ]; then
	grep ">" wide/MseI_MspI.fasta | \
	gawk '{split($1,A,/[-\|]/); F[A[4]-A[3]]++}END{for (f in F) print f "\t" F[f]}' | \
	sort -nk 1 > wide/lengths_all.hist
fi

if [ ! -e plot.gnp ]; then
	echo "set output 'lengths.png'" > plot.gnp
	echo "set term png large" >> plot.gnp
	echo "set ylab 'Frequency'" >> plot.gnp
	echo "set xlab 'Length (bp)'" >> plot.gnp
	echo "plot 'wide/lengths_all.hist' using 1:2 with lines title 'all fragments', 'wide/lengths_inner.hist' using 1:2 with lines title 'inner fragments'" >> plot.gnp
	echo "quit" >> plot.gnp
fi

if [ ! -e lengths.png ]; then
	gnuplot < plot.gnp
fi

# Number of fragments in every 50 bp size range between 450 and 749.
if [ ! -e wide/fragNumRange50_all.txt ]; then
	gawk '{F[$1] = $2}END{for (i=450; i<=700; i++) {S = 0; for (j=i; j<i+50; j++) S += F[j]; print i "\t" j "\t" S}}' wide/lengths_all.hist > wide/fragNumRange50_all.txt
fi

if [ ! -e wide/fragNumRange50_inner.txt ]; then
        gawk '{F[$1] = $2}END{for (i=450; i<=700; i++) {S = 0; for (j=i; j<i+50; j++) S += F[j]; print i "\t" j "\t" S}}' wide/lengths_inner.hist > wide/fragNumRange50_inner.txt
fi

if [ ! -e plot2.gnp ]; then
	echo "set output 'fragments.png'" > plot2.gnp
	echo "set term png large" >> plot2.gnp
	echo "set xlab 'Minimum size (bp)'" >> plot2.gnp
	echo "set ylab 'Number of fragments'" >> plot2.gnp
	echo "set title '50-bp ranges'" >> plot2.gnp
	echo "plot 'wide/fragNumRange50_all.txt' using 1:3 with lines title 'all fragments', 'wide/fragNumRange50_inner.txt' using 1:3 with lines title 'inner fragments'" >> plot2.gnp
	echo "quit" >> plot2.gnp
fi

if [ ! -e fragments.png ]; then
	gnuplot < plot2.gnp
fi

# Results
#
# The plot lenghts.png shows the distribution of fragment lengths including (all fragments) or
# excluding (inner fragments) the fragments flanked by contig ends. The size range of interest
# is flanked by two peaks, at 494 and 608 bp. I believe these two peaks to be homologous to
# the peaks observed at around 450 and 510 in Coregonus samples (table below). Although the
# pairs of peaks do not perfectly coincide between the two species, they cover a similar range
# and have similar relative magnitudes: the first peak is much higher than the second. In fact,
# the avoidance of the second peak is not much of a concern in comparison to the first one.
#
# ----------------------------
# Sample	peak	peak
# 		(bp)	(bp)
# ---------------------------
# St001		446	514
# St003		445	511
# St006		446	512
# St015		450	516
# St016		453	520
# St019		454	516
# St037		453	523
# St039		454	518
# St043		453	518
# St049		453	517
# St050		452	514
# Bl065		445	517
# Bl076		444	510
# Bl080		442	510
# Bl083		447	510
# Bl104		447	511
# Bl108		447	507
# Bl091		448	511
# Bl094		447	509
# Bl095		450	505
# Bl098		446	508
# Bl116		447	510
# ---------------------------------------
#
# I determine the number of fragments to be slightly high, in the order of 40000 to
# 80000. That means that with 17 million reads, and 12 samples per run, we would get
# between 17 and 35 fold coverage per site and individual. Not very large for
# genotyping.
