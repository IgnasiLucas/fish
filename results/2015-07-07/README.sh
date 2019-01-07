#!/bin/bash
#
#				2015-07-07
#				----------
#
# The first batch of ligation reactions gave results that suggest extensive
# re-ligation of genomic fragments. I suspect low quality of adapters. I will
# re-anneal adapters, check their concentration, and repeat the ligations in
# samples that have enough DNA left over to attempt another ligation.

DATADIR=`pwd | sed "s/results/data/"`
if [ ! -d $DATADIR ]; then
	mkdir $DATADIR
fi

# This time, I will run only one reaction per sample.

Sample=(0 St001 St003 St006 St015 St016 St019 St037 St039 St043 St044 St049 St050 \
	  Bl065 Bl076 Bl080 Bl083 Bl104 Bl108 Bl091 Bl093 Bl094 Bl095 Bl098 Bl116)

# These are the volumes of sample DNA that I will use in the second ligation. They make
# up to 500 ng per sample, although several samples do not reach that amount, and some
# do not have any DNA left.

Volume=(0 32.8 43.0  0.0 43.2 50.1 45.7 44.0 50.8 49.2 43.1 61.0 48.9 \
          20.3 33.5 53.7  0.0  0.0 49.2 59.9 55.1 53.2 57.8  0.0 40.9)

# I will use 8 P5 adapters, with different codewords. The P7 adapter is common to all samples.
# Thus, I will only be able to make three 8-samples pools. Each pools will be kept separate,
# and it will get one of three indexes in the P7 side after the amplification PCR. The easy
# way to distribute the 8 P5 adapters is sequentially. However, it is desirable to balance
# the pools and sequencing runs across the 4 sample groups. Thus, each pool of 8 samples should
# have two individuals from each group. For example:
#
#      ---------------------------------------------------------------------
#       Sample  	Group                               Adapter    Pool
#      ---------------------------------------------------------------------
# 	St0001  	Coregonus albula/Stechlinsee    	A01	1
# 	St0003  	Coregonus albula/Stechlinsee    	B01	1
# 	St0006  	Coregonus albula/Stechlinsee    	G01	2
# 	St0015  	Coregonus albula/Stechlinsee    	H01	2
# 	St0016-2	Coregonus albula/Stechlinsee    	E01	3
# 	St0019-2	Coregonus albula/Stechlinsee    	F01	3
# 	St0037-2	Coregonus fontanae/Stechlinsee   	C01	1
# 	St0039-2	Coregonus fontanae/Stechlinsee   	D01	1
# 	St0043-2	Coregonus fontanae/Stechlinsee    	A01	2
# 	St0044-2	Coregonus fontanae/Stechlinsee    	B01	2
# 	St0049-2	Coregonus fontanae/Stechlinsee    	G01	3
# 	St0050-2	Coregonus fontanae/Stechlinsee    	H01	3
# 	Bl0065-2	Coregonus albula/Breiter Luzin    	E01	1
# 	Bl0076  	Coregonus albula/Breiter Luzin    	F01	1
# 	Bl0080  	Coregonus albula/Breiter Luzin    	C01	2
# 	Bl0083-2	Coregonus albula/Breiter Luzin    	D01	2
# 	Bl0104  	Coregonus albula/Breiter Luzin    	A01	3
# 	Bl0108-2	Coregonus albula/Breiter Luzin    	B01	3
# 	Bl0091-2	Coregonus lucinensis/Breiter Luzin	G01	1
# 	Bl0093-2	Coregonus lucinensis/Breiter Luzin	H01	1
# 	Bl0094-2	Coregonus lucinensis/Breiter Luzin	E01	2
# 	Bl0095-2	Coregonus lucinensis/Breiter Luzin	F01	2
# 	Bl0098-2	Coregonus lucinensis/Breiter Luzin	C01	3
# 	Bl0116-2	Coregonus lucinensis/Breiter Luzin	D01	3
#      -------------------------------------------------------------------
#

Adapter=(0 A01 B01 G01 H01 E01 F01 C01 D01 A01 B01 G01 H01 \
           E01 F01 C01 D01 A01 B01 G01 H01 E01 F01 C01 D01)

Pool=(0 1 1 2 2 3 3 1 1 2 2 3 3 \
        1 1 2 2 3 3 1 1 2 2 3 3)

# The following are average fragment masses, in ng/fmol. They can be calculated from the Bioanalyzer
# results in two ways. First, if you divide concentration (pg/µl) by molarity (pmol/l, which
# is equivalent to attomol/µl), you get the average fragment mass in pg/amol, which is
# equivalent to ng/fmol. Otherwise, you can use the reported average fragment size and apply
# the formula of the molecular mass of double-stranded DNA, namely  (# nucleotides x 607.4) + 157.9,
# to obtain the g/mol, and divide by 10^6 to express it in ng/fmol. However, the results are quite
# different. The second method gives a 3-fold higher value. The ligation molarity calculator
# cited by the ddRAD protocol uses the first method. I don't understand the difference, but I
# will just use the method recommended by the ligation molarity calculator. Find the data in
# doc/samplelist_20150420.ods. I assign the average fragment mass across the rest of the samples
# to samples without BioAnalyzer data (St044 and Bl093), namely 0.210487 ng/fmol.

FragMass=(0 0.218318 0.206906 0.202655 0.207858 0.212878 0.205778 0.210686 0.213317 0.216652 0.211076 0.218938 0.214750 \
            0.219916 0.211229 0.233630 0.202335 0.210486 0.207836 0.217333 0.209892 0.210238 0.204225 0.211104 0.190475)

# These are the DNA concentrations (ng/µl), obtained from Bioanalyzer.
Concentration=(0 9.5 11.6 10.9 11.6 10.0 10.9 11.4  9.8 10.2 10.4  8.2 10.2 \
                 8.7  8.7  9.3  7.6 10.3 10.2  8.4  9.1  9.4  8.7  8.5 12.2)

# The original adapter concentration is assumed to be equal for both adapters: 40 ng/µl.

StockConc=40

# The DNA is fragmented with two enzymes. There are three types of fragments: MseI-MseI, MspI-MspI,
# and either MseI-MspI or MspI-MseI. The fact is that I do not know the proportion of each type
# of fragment, and therefore I cannot calculate the molarity of each type of fragment end very
# accurately. One option is to calculate the maximum molarity of each type, as if all ends were
# of the same type. Another option is to use the proportion of fragment ends observed in the in
# silico digestion of Salmo salar's reference genome (2015-06-15), which is:
#
# MseI (TTAA): 80.58%	(AT content: )
# MspI (CCGG): 19.42%	(CG content: )
#
# In the first case, I would set the fold-excess to no more than 7, since I would already be
# overestimating the molarity of fragment ends. In the second case, I would set the fold-excess of
# adapters to ends to 10. Given how unbalanced the two types of ends are, I think it makes sense
# to account for it on the bases of Salmo salar's genome.
#
# Recall: MseI produces 5'-TA overhangs, that bind to barcoded P5 adapter.
#         MspI produces 5'-CG overhangs, that bind to common, unlabeled P7 adapter.
#
# The original ddRAD protocol prepares two different adapter working stocks, one for each
# adapter. Here, I am trying to prepare one solution with both adapters, each at its right
# concentration. Is it too ambitious? The added difficulty is that the dilution factor of each
# adapter from the stock to the working solution must be large enough to accommodate the required
# volume of the other adapter.

echo -e "Sample\tP5 adapter\tPool\tDNA\tAdapters\t10X T4 buffer\t1.5M NaCl\tT4 Ligase\tWater\tTotal" > ligations.txt
echo -e "Sample\tDNA\tP5 Adapter\tP5 (µl)\tP7 (µl)\tAnnealing buffer (µl)\tTotal" > workingStock.txt
if [ -e results ]; then
	rm results
fi
for i in `seq 1 24`; do
        perl -we 'use POSIX;  # "use POSIX" allows me to use the "ceil" function below.
                my $Conc = $ARGV[0]; # DNA concentration, ng/µl.
		my $Vol = $ARGV[1];  # DNA volume, µl.
		my $FragMass = $ARGV[2]; # Average DNA fragment mass, ng/fmol.
		my $PropP5 = 0.8058; # Expected proportion of fragment ends that stick to P5 adapters.
		my $PropP7 = 0.1942; # Expected proportion of fragment ends that stick to P7 adapters.
                my $AdapterStockConc = $ARGV[3]; # Annealed adapter stocks concentration (µM, or pmol / µl). Assumed equal for both adapters.
                my $Ends = 2 * $Conc * $Vol / ($FragMass * 1000); # Total number of fragment ends per sample (pmol).
                my $FoldExcess = 12;
                my $Adapter_P5 = $Ends * $PropP5 * $FoldExcess; # Amount of P5 adapter per sample (pmols).
		my $Adapter_P7 = $Ends * $PropP7 * $FoldExcess; # Amount of P7 adapter per sample (pmols).
                my $AdapterVolRxn = 1; # Target adapter volum/reaction (µl).
                my $AdapterWorkConc_P5 = $Adapter_P5 / $AdapterVolRxn; # Required concentration of adapter P5 in working stock (pmols / µl or µM).
		my $AdapterWorkConc_P7 = $Adapter_P7 / $AdapterVolRxn; # Required concentration of adpater P7 in working stock (pmols / µl or µM).
		# Below I require the dilution factor of each adapter to be at least equal to the expected
		# proportion of the corresponding type of fragment ends. Because the two adapters will be
		# mixed in the same dilution, and each one will become diluted at least to its mixing proportion.
		while (($AdapterWorkConc_P5 > $PropP5 * $AdapterStockConc) || ($AdapterWorkConc_P7 > $PropP7 * $AdapterStockConc)) {
			$AdapterVolRxn += 0.25;
			$AdapterWorkConc_P5 = $Adapter_P5 / $AdapterVolRxn;
			$AdapterWorkConc_P7 = $Adapter_P7 / $AdapterVolRxn;
		}
                my $AdapterStockVol_P5 = 0; # Volume of adapter P5 stock (at $AdapterStockConc) needed to prepare $AdapterWorkVol µl at $AdapterWorkConc_P5 (µl).
		my $AdapterStockVol_P7 = 0; # Volume of adapter P7 stock (at $AdapterStockConc) needed to prepare $AdapterWorkVol µl at $AdapterWorkConc_P7 (µl).
						# $AdapterStockVol_P5 * $AdapterStockConc == $AdapterWorkVol * $AdapterWorkConc_P5
						# $AdapterStockVol_P7 * $AdapterStockConc == $AdapterWorkVol * $AdapterWorkConc_P7
                my $AdapterWorkVol = 7; # Working stock volum of annealed adapters to make (µl).
                # Below, I make sure that I do not have to pipette less than 0.5 µl of each adapter stock.
		if ($Vol > 0) {
			while (($AdapterStockVol_P5 < 1) || ($AdapterStockVol_P7 < 1) || ($AdapterStockVol_P5 + $AdapterStockVol_P7 > $AdapterWorkVol)) {
	                        $AdapterWorkVol++;
	                        $AdapterStockVol_P5 = $AdapterWorkVol * $AdapterWorkConc_P5 / $AdapterStockConc;
				$AdapterStockVol_P7 = $AdapterWorkVol * $AdapterWorkConc_P7 / $AdapterStockConc;
	                }
		} else {
			$AdapterWorkVol = 0;
			$AdapterVolRxn = 0;
		}
		my $AnnealBufferVol = $AdapterWorkVol - $AdapterStockVol_P5 - $AdapterStockVol_P7;
		# The reaction includes:
		# 	$Vol µl of DNA,
		#	$AdapterVolRxn µl of adapters,
		#	(1/10) * $TotVol µl of 10X T4 ligase buffer,
		#	(1/40) * $TotVol µl of 1.5M NaCl (which is ~40X, relative to the intended final concentration of 38 mM),
		#	1 µl of T4 DNA ligase,
		#	$TotVol - (all the above) µl of water.
                my $TotVol = ceil(($Vol + 1 + $AdapterVolRxn) / (1 - 1/10 - 1/40));
                print "Sample $ARGV[4]:\n";
                print "  Amount of DNA: $Vol µl at $Conc ng/µl.\n";
                print "  Average fragment mass (from BioAnalyzer): $FragMass ng/fmol.\n";
                print "  Fragment ends: ", $Ends, " pmols.\n";
		print "     Expected amount of P5-compatible ends: ", $Ends * $PropP5, " pmols.\n";
		print "     Expected amount of P7-compatible ends: ", $Ends * $PropP7, " pmols.\n";
                print "  Target amount of P5 adapter $ARGV[5] (P5-compatible ends \* $FoldExcess): ", $Adapter_P5, " pmols.\n";
		print "  Target amount of common adapter P7 (P7-compatible ends \* $FoldExcess): ", $Adapter_P7, " pmols.\n";
                print "  Concentration of P5 adapter $ARGV[5] required: ", $AdapterWorkConc_P5, " pmol/µl.\n";
		print "  Concentration of common P7 adapter required: ", $AdapterWorkConc_P7, " pmol/µl.\n";
                print "  Concentration of common and $ARGV[5] adapter stocks: ", $AdapterStockConc, " pmol/µl.\n";
                print "  Volume of working adapter mixture stock to make: ", $AdapterWorkVol, " µl.\n";
                print "  Volume of original P5 adapter $ARGV[5] stock required: ", $AdapterStockVol_P5, " µl.\n";
		print "  Volume of original P7 common adapter stock required: ", $AdapterStockVol_P7, " µl.\n";
                print "  Volume of 1X annealing buffer: ", $AnnealBufferVol, " µl.\n";
                print "  Volume of this working stock per reaction: ", $AdapterVolRxn, " µl.\n";
                print "  Ligation reaction $ARGV[4]: $Vol µl DNA, $AdapterVolRxn µl adapter mixture working stock, ";
                print sprintf("%.2f", $TotVol/10), " µl 10X Ligase Buffer, ", sprintf("%.2f", $TotVol/40), " µl 1.5M NaCl, ";
                print "1 µl T4 Ligase, and water to $TotVol µl.\n\n";

		open(LIG, ">>", "ligations.txt") || die "I cannot open ligations.txt.\n";
		open(WRK, ">>", "workingStock.txt") || die "I cannot open workingSotck.txt.\n";
		my $water = $TotVol - $Vol - $AdapterVolRxn - $TotVol/10 - $TotVol/40 - 1;
		print LIG $ARGV[4], "\t", $ARGV[5], "\t", $ARGV[6], "\t", $Vol, "\t", $AdapterVolRxn, "\t", sprintf("%.2f", $TotVol/10), "\t", sprintf("%.2f", $TotVol/40), "\t1.0\t", sprintf("%.2f", $water), "\t", $TotVol, "\n";
		print WRK $ARGV[4], "\t", $ARGV[5], "\t", sprintf("%.2f\t%.2f\t%.2f\t%.2f\n", $AdapterStockVol_P5, $AdapterStockVol_P7, $AnnealBufferVol, $AdapterWorkVol);
		close(LIG);
		close(WRK);' ${Concentration[$i]} ${Volume[$i]} ${FragMass[$i]} $StockConc ${Sample[$i]} ${Adapter[$i]} ${Pool[$i]} >> results
done

# In the calculations above, I set the volume of adapter working stock as low as possible, and
# determined the composition of the working stock according to that volume and to the sample's
# molarity. As a result, a different working stock is prescribed for each sample, even among
# samples sharing the same P5 adapter. In fact, all the working stocks have the same proportions
# of the two adapters. Thus, it must be possible to prepare only 8 working stocks (one per P5
# adapter), instead of 24 (one per sample). I can fix the volume of adapter working stock to
# a reasonable value (7 µl?), and allow the volume of adapter working stock per reaction to
# vary.

echo -e "Sample\tP5 adapter\tPool\tDNA\tAdapters\t10X T4 buffer\t1.5M NaCl\tT4 Ligase\tWater\tTotal" > ligations2.txt
echo -e "Sample\tP5 Adapter\tP5 (µl)\tP7 (µl)\tAnnealing buffer (µl)\tTotal (µl)" > workingStock2.txt
if [ -e results_2 ]; then
	rm results_2
fi
for i in `seq 1 24`; do
        perl -we 'use POSIX;  # "use POSIX" allows me to use the "floor" function below.
                my $Conc = $ARGV[0]; # DNA concentration, ng/µl.
                my $Vol = $ARGV[1];  # DNA volume, µl.
                my $FragMass = $ARGV[2]; # Average DNA fragment mass, ng/fmol.
                my $PropP5 = 0.8; #0.8058; # Expected proportion of fragment ends that stick to P5 adapters.
                my $PropP7 = 0.2; #0.1942; # Expected proportion of fragment ends that stick to P7 adapters.
                my $AdapterStockConc = $ARGV[3]; # Annealed adapter stocks molarity (µM, or pmol / µl). Assumed equal for both adapters.
                my $Ends = 2 * $Conc * $Vol / ($FragMass * 1000); # Total number of fragment ends per sample (pmol).
                my $FoldExcess = 12;
                my $Adapter_P5 = $Ends * $PropP5 * $FoldExcess; # Amount of P5 adapter per sample (pmols).
                my $Adapter_P7 = $Ends * $PropP7 * $FoldExcess; # Amount of P7 adapter per sample (pmols).
		my $AdapterWorkVol = 7; # Predetermined volume of adapter working stock, to re-use working stocks across pools (µl).
		# It is easy to realize that $AdapterStockVol_P5 / $AdapterStockVol_P7 = Adapter_P5 / Adapter_P7.
		# In addition, $AdapterStockVol_P5 + $AdapterStockVol_P7 <= $AdapterWorkVol, the rest being $AnnealBufferVol.
		# It follows that $AdapterStockVol_P7 <= $AdapterWorkVol / (1 + PropP5 / PropP7). Thus, I set the following:
		my $AdapterStockVol_P7 = floor($AdapterWorkVol / (1 + $PropP5 / $PropP7)); # Volume of common P7 adapter stock used to make the working stock (µl).
		my $AdapterStockVol_P5 = $AdapterStockVol_P7 * $PropP5 / $PropP7;	   # Volume of labeled P5 adapter stock used to make the working stock (µl).
		my $AnnealBufferVol = $AdapterWorkVol - $AdapterStockVol_P5 - $AdapterStockVol_P7; # Volume of annealing buffer required to make the working stock (µl).
		my $AdapterWorkConc_P5 = $AdapterStockVol_P5 * $AdapterStockConc / $AdapterWorkVol; # Molarity of P5 in working stock (pmol/µl).
		my $AdapterWorkConc_P7 = $AdapterStockVol_P7 * $AdapterStockConc / $AdapterWorkVol; # Molarity of P7 in working stock (pmol/µl)
		my $AdapterVolRxn = $Adapter_P5 / $AdapterWorkConc_P5;
		if (sprintf("%.4f", $AdapterVolRxn) ne sprintf("%.4f", $Adapter_P7 / $AdapterWorkConc_P7)) {
			print "Warning: The adapter volume per reactions is different if calculated\n";
			print "         from the other adapter:\n";
			print "         Adapter_P7 / AdapterWorkConc_P7 = ", $Adapter_P7 / $AdapterWorkConc_P7, "\n";
			print "         Adapter_P5 / AdapterWorkConc_P5 = ", $Adapter_P5 / $AdapterWorkConc_P5, " (using this).\n";
		}

                my $TotVol = ceil(($Vol + 1 + $AdapterVolRxn) / (1 - 1/10 - 1/40));
                print "Sample $ARGV[4]:\n";
                print "  Amount of DNA: $Vol µl at $Conc ng/µl.\n";
                print "  Average fragment mass (from BioAnalyzer): $FragMass ng/fmol.\n";
                print "  Fragment ends: ", $Ends, " pmols.\n";
                print "     Expected amount of P5-compatible ends: ", $Ends * $PropP5, " pmols.\n";
                print "     Expected amount of P7-compatible ends: ", $Ends * $PropP7, " pmols.\n";
                print "  Target amount of P5 adapter $ARGV[5] (P5-compatible ends \* $FoldExcess): ", $Adapter_P5, " pmols.\n";
                print "  Target amount of common adapter P7 (P7-compatible ends \* $FoldExcess): ", $Adapter_P7, " pmols.\n";
                print "  Concentration of P5 adapter $ARGV[5] required: ", $AdapterWorkConc_P5, " pmol/µl.\n";
                print "  Concentration of common P7 adapter required: ", $AdapterWorkConc_P7, " pmol/µl.\n";
                print "  Concentration of common and $ARGV[5] adapter stocks: ", $AdapterStockConc, " pmol/µl.\n";
                print "  Volume of working adapter mixture stock to make: ", $AdapterWorkVol, " µl.\n";
                print "  Volume of original P5 adapter $ARGV[5] stock required: ", $AdapterStockVol_P5, " µl.\n";
                print "  Volume of original P7 common adapter stock required: ", $AdapterStockVol_P7, " µl.\n";
                print "  Volume of 1X annealing buffer: ", $AnnealBufferVol, " µl.\n";
                print "  Volume of this working stock per reaction: ", $AdapterVolRxn, " µl.\n";
                print "  Ligation reaction $ARGV[4]: $Vol µl DNA, $AdapterVolRxn µl adapter mixture working stock, ";
                print sprintf("%.2f", $TotVol/10), " µl 10X Ligase Buffer, ", sprintf("%.2f", $TotVol/40), " µl 1.5M NaCl, ";
                print "1 µl T4 Ligase, and water to $TotVol µl.\n\n";

                open(LIG, ">>", "ligations2.txt") || die "I cannot open ligations2.txt.\n";
                open(WRK, ">>", "workingStock2.txt") || die "I cannot open workingSotck2.txt.\n";
                my $water = $TotVol - $Vol - $AdapterVolRxn - $TotVol/10 - $TotVol/40 - 1;
                print LIG $ARGV[4], "\t", $ARGV[5], "\t", $ARGV[6], "\t", $Vol, "\t", sprintf("%.2f\t%.2f\t%.2f\t1.0\t%.2f\t%.2f\n", $AdapterVolRxn, $TotVol/10, $TotVol/40, $water, $TotVol);
                print WRK $ARGV[4], "\t", $ARGV[5], "\t", sprintf("%.2f\t%.2f\t%.2f\t%.2f\n", $AdapterStockVol_P5, $AdapterStockVol_P7, $AnnealBufferVol, $AdapterWorkVol);
                close(LIG);
                close(WRK);' ${Concentration[$i]} ${Volume[$i]} ${FragMass[$i]} $StockConc ${Sample[$i]} ${Adapter[$i]} ${Pool[$i]} >> results_2
done
