#!/bin/bash
#
#				2015-05-05
#				----------
#
# This is to calculate the composition of the ligation reactions. In this project, I am
# using two adapters, the P5 adapter is labeled, while the P7 adapter is common to all
# samples. I want to prepare an adapter working stock with both adapters in it, each at
# the desired concentration. This is only possible if I prepare a volume of working stock
# that is at least as much as the joint volumes of the two original adapter stocks
# required. In turn, this is only feasible if the desired concentration of the working
# stock is not too high.

DATADIR=`pwd | sed "s/results/data/"`
if [ ! -d $DATADIR ]; then
	mkdir $DATADIR
fi

# The following is some data about the samples from doc/samplelist_20150420.ods and
# from the Bioanalyzer results.
if [ ! -e $DATADIR/data ]; then
	echo -e "\tSp\tLake\tMass\tNum.Digest\tVol.digest\tVol.total\tElution\tAve.Size\tConc\tMolarity" > $DATADIR/data
	echo -e "St0001\tCoregonus albula\tStechlinsee\t2910.6\t3\t33.0\t132.00\t97.02\t1082\t3162.08\t14483.8" >> $DATADIR/data
	echo -e "St0003\tCoregonus albula\tStechlinsee\t7177.5\t7\t13.7\t54.69\t239.25\t1058\t3874.93\t18728.0" >> $DATADIR/data
	echo -e "St0006\tCoregonus albula\tStechlinsee\t1881.0\t2\t49.5\t198.00\t62.70\t1048\t3628.22\t17903.4" >> $DATADIR/data
	echo -e "St0015\tCoregonus albula\tStechlinsee\t5329.0\t5\t29.3\t117.12\t177.63\t1128\t3860.07\t18570.7" >> $DATADIR/data
	echo -e "St0016-2\tCoregonus albula\tStechlinsee\t3801.6\t4\t24.8\t99.00\t126.72\t1213\t3328.79\t15637.1" >> $DATADIR/data
	echo -e "St0019-2\tCoregonus albula\tStechlinsee\t7929.0\t8\t11.9\t47.65\t264.30\t1071\t3646.72\t17721.6" >> $DATADIR/data
	echo -e "St0037-2\tCoregonus fontanae\tStechlinsee\t5542.7\t6\t16.2\t64.60\t184.76\t1137\t3788.53\t17981.9" >> $DATADIR/data
	echo -e "St0039-2\tCoregonus fontanae\tStechlinsee\t4772.5\t5\t19.6\t78.56\t159.08\t1148\t3278.98\t15371.4" >> $DATADIR/data
	echo -e "St0043-2\tCoregonus fontanae\tStechlinsee\t3643.2\t4\t24.8\t99.00\t121.44\t1171\t3384.14\t15620.2" >> $DATADIR/data
	echo -e "St0044-2\tCoregonus fontanae\tStechlinsee\t3029.4\t3\t33.0\t132.00\t100.98\tNA\tNA\tNA" >> $DATADIR/data
	echo -e "St0049-2\tCoregonus fontanae\tStechlinsee\t5633.8\t6\t16.1\t64.53\t187.79\t1223\t2731.03\t12474.0" >> $DATADIR/data
	echo -e "St0050-2\tCoregonus fontanae\tStechlinsee\t3286.8\t3\t33.0\t132.00\t109.56\t1169\t3404.89\t15855.1" >> $DATADIR/data
	echo -e "Bl0065-2\tCoregonus albula\tBreiter Luzin\t2712.6\t3\t33.0\t132.00\t90.42\t1233\t2893.90\t13159.1" >> $DATADIR/data
	echo -e "Bl0076\tCoregonus albula\tBreiter Luzin\t3099.2\t3\t49.7\t198.67\t103.31\t1143\t2905.48\t13755.1" >> $DATADIR/data
	echo -e "Bl0080\tCoregonus albula\tBreiter Luzin\t4052.8\t4\t37.3\t149.00\t135.09\t1241\t3105.72\t13293.3" >> $DATADIR/data
	echo -e "Bl0083-2\tCoregonus albula\tBreiter Luzin\t1098.9\t2\t49.5\t198.00\t36.63\t1160\t2536.77\t12537.5" >> $DATADIR/data
	echo -e "Bl0104\tCoregonus albula\tBreiter Luzin\t1910.7\t2\t49.5\t198.00\t63.69\t1139\t3421.39\t16254.7" >> $DATADIR/data
	echo -e "Bl0108-2\tCoregonus albula\tBreiter Luzin\t8024.3\t8\t11.9\t47.65\t267.48\t1144\t3390.85\t16315.0" >> $DATADIR/data
	echo -e "Bl0091-2\tCoregonus lucinensis\tBreiter Luzin\t5782.7\t6\t16.1\t64.47\t192.76\t1237\t2783.73\t12808.6" >> $DATADIR/data
	echo -e "Bl0093-2\tCoregonus lucinensis\tBreiter Luzin\t3841.2\t4\t24.8\t99.00\t128.04\tNA\tNA\tNA" >> $DATADIR/data
	echo -e "Bl0094-2\tCoregonus lucinensis\tBreiter Luzin\t4865.8\t5\t19.6\t78.48\t162.19\t1180\t3135.47\t14913.9" >> $DATADIR/data
	echo -e "Bl0095-2\tCoregonus lucinensis\tBreiter Luzin\t4197.6\t4\t24.8\t99.00\t139.92\t1158\t2885.56\t14129.3" >> $DATADIR/data
	echo -e "Bl0098-2\tCoregonus lucinensis\tBreiter Luzin\t2118.6\t2\t49.5\t198.00\t70.62\t1153\t2816.78\t13343.1" >> $DATADIR/data
	echo -e "Bl0116-2\tCoregonus lucinensis\tBreiter Luzin\t3484.8\t3\t33.0\t132.00\t116.16\t1064\t4073.27\t21384.8" >> $DATADIR/data
fi

R --save < plot.R

# I will split each sample in two reactions, to keep volumes low. However, sample Bl083
# did not get enough DNA, and I have enough with only one reaction. Thus, there must be
# 24 * 2 - 1 = 47 reactions.

Sample=(0 St001.1 St001.2 St003.1 St003.2 St006.1 St006.2 St015.1 St015.2 St016.1 St016.2 St019.1 St019.2 \
	  St037.1 St037.2 St039.1 St039.2 St043.1 St043.2 St044.1 St044.2 St049.1 St049.2 St050.1 St050.2 \
	  Bl065.1 Bl065.2 Bl076.1 Bl076.2 Bl080.1 Bl080.2 Bl083.1         Bl104.1 Bl104.2 Bl108.1 Bl108.2 \
	  Bl091.1 Bl091.2 Bl093.1 Bl093.2 Bl094.1 Bl094.2 Bl095.1 Bl095.2 Bl098.1 Bl098.2 Bl116.1 Bl116.2)

# These are the volumes of sample DNA required to use up to 600 ng of each sample (except
# sample Bl083 and Bl098, which have 270 and 588 ng, respectively). See spreadsheet for
# details (doc/samplelist_20150420.ods).

Volume=(0 31.62 31.62 25.81 25.81 27.56 27.56 25.91 25.91 30.04 30.04 27.42 27.42 \
	  26.40 26.40 30.50 30.50 29.55 29.55 30.00 30.00 36.62 36.62 29.37 29.37 \
	  34.56 34.56 34.42 34.42 32.20 32.20 35.64       29.23 29.23 29.49 29.49 \
	  35.92 35.92 37.50 37.50 31.89 31.89 34.66 34.66 34.81 34.81 24.55 24.55)

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
# Now, keeping in mind that each sample is divided in two aliquotes, except
# Bl083, then:

Adapter=(0 A01 A01 B01 B01 G01 G01 H01 H01 E01 E01 F01 F01 \
	   C01 C01 D01 D01 A01 A01 B01 B01 G01 G01 H01 H01 \
	   E01 E01 F01 F01 C01 C01 D01     A01 A01 B01 B01 \
	   G01 G01 H01 H01 E01 E01 F01 F01 C01 C01 D01 D01)

Pool=(0 1 1 1 1 2 2 2 2 3 3 3 3
        1 1 1 1 2 2 2 2 3 3 3 3
        1 1 1 1 2 2 2   3 3 3 3
        1 1 1 1 2 2 2 2 3 3 3 3)

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

FragMass=(0 0.218318 0.218318 0.206906 0.206906 0.202655 0.202655 0.207858 0.207858 0.212878 0.212878 0.205778 0.205778 \
	    0.210686 0.210686 0.213317 0.213317 0.216652 0.216652 0.210487 0.210487 0.218938 0.218938 0.214750 0.214750 \
	    0.219916 0.219916 0.211229 0.211229 0.233630 0.233630 0.202335          0.210486 0.210486 0.207836 0.207836 \
	    0.217333 0.217333 0.210487 0.210487 0.210238 0.210238 0.204225 0.204225 0.211104 0.211104 0.190475 0.190475)

# These are the DNA concentrations (ng/µl), obtained from Bioanalyzer.

Concentration=(0 9.5  9.5 11.6 11.6 10.9 10.9 11.6 11.6 10.0 10.0 10.9 10.9 \
		11.4 11.4  9.8  9.8 10.2 10.2 10.4 10.4  8.2  8.2 10.2 10.2 \
		 8.7  8.7  8.7  8.7  9.3  9.3  7.6      10.3 10.3 10.2 10.2 \
		 8.4  8.4  9.1  9.1  9.4  9.4  8.7  8.7  8.5  8.5 12.2 12.2)

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
for i in `seq 1 47`; do
        perl -we 'use POSIX;  # "use POSIX" allows me to use the "ceil" function below.
                my $Conc = $ARGV[0]; # DNA concentration, ng/µl.
		my $Vol = $ARGV[1];  # DNA volume, µl.
		my $FragMass = $ARGV[2]; # Average DNA fragment mass, ng/fmol.
		my $PropP5 = 0.8058; # Expected proportion of fragment ends that stick to P5 adapters.
		my $PropP7 = 0.1942; # Expected proportion of fragment ends that stick to P7 adapters.
                my $AdapterStockConc = $ARGV[3]; # Annealed adapter stocks concentration (µM, or pmol / µl). Assumed equal for both adapters.
                my $Ends = 2 * $Conc * $Vol / ($FragMass * 1000); # Total number of fragment ends per sample (pmol).
                my $FoldExcess = 10;
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
		while (($AdapterStockVol_P5 < 2) || ($AdapterStockVol_P7 < 2) || ($AdapterStockVol_P5 + $AdapterStockVol_P7 > $AdapterWorkVol)) {
                        $AdapterWorkVol++;
                        $AdapterStockVol_P5 = $AdapterWorkVol * $AdapterWorkConc_P5 / $AdapterStockConc;
			$AdapterStockVol_P7 = $AdapterWorkVol * $AdapterWorkConc_P7 / $AdapterStockConc;
                }
                my $AnnealBufferVol = $AdapterWorkVol - $AdapterStockVol_P5 - $AdapterStockVol_P7;
		# The reaction includes:
		# 	$Vol µl of DNA,
		#	$AdapterVolRxn µl of adapters,
		#	(1/10) * $TotVol µl of 10X T4 ligase buffer,
		#	(1/40) * $TotVol µl of 1.5M NaCl (which is ~40X, relative to the intended final concentration of 38 mM),
		#	2 µl of T4 DNA ligase,
		#	$TotVol - (all the above) µl of water.
                my $TotVol = ceil(($Vol + 2 + $AdapterVolRxn) / (1 - 1/10 - 1/40));
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
                print "2 µl T4 Ligase, and water to $TotVol µl.\n\n";

		open(LIG, ">>", "ligations.txt") || die "I cannot open ligations.txt.\n";
		open(WRK, ">>", "workingStock.txt") || die "I cannot open workingSotck.txt.\n";
		my $water = $TotVol - $Vol - $AdapterVolRxn - $TotVol/10 - $TotVol/40 - 2;
		print LIG $ARGV[4], "\t", $ARGV[5], "\t", $ARGV[6], "\t", $Vol, "\t", $AdapterVolRxn, "\t", sprintf("%.2f", $TotVol/10), "\t", sprintf("%.2f", $TotVol/40), "\t2.0\t", sprintf("%.2f", $water), "\t", $TotVol, "\n";
		print WRK $ARGV[4], "\t", $ARGV[5], "\t", sprintf("%.2f\t%.2f\t%.2f\t%.2f\n", $AdapterStockVol_P5, $AdapterStockVol_P7, $AnnealBufferVol, $AdapterWorkVol);
		close(LIG);
		close(WRK);' ${Concentration[$i]} ${Volume[$i]} ${FragMass[$i]} $StockConc ${Sample[$i]} ${Adapter[$i]} ${Pool[$i]} >> results
done

# In the calculations above, I set the volume of adapter working stock as low as possible, and
# determined the composition of the working stock according to that volume and to the sample's
# molarity. As a result, a different working stock is prescribed for each sample, even among
# samples sharing the same P5 adapter. In fact, all the working stocks have the same proportions
# of the two adapters. Thus, it must be possible to prepare only 8 working stocks (one per P5
# adapter), instead of 24 (one per sample). In view of the results from above, I know that the
# largest amount of adapter working stock is 17 µl. Now , I should be able to re-work the
# calculations fixing all volumes of adapter working stocks to 17 µl, and allowing for variation
# in the volume of adapter working stock per reaction.
#

echo -e "Sample\tP5 adapter\tPool\tDNA\tAdapters\t10X T4 buffer\t1.5M NaCl\tT4 Ligase\tWater\tTotal" > ligations2.txt
echo -e "Sample\tP5 Adapter\tP5 (µl)\tP7 (µl)\tAnnealing buffer (µl)\tTotal (µl)" > workingStock2.txt
if [ -e results_2 ]; then
	rm results_2
fi
for i in `seq 1 47`; do
        perl -we 'use POSIX;  # "use POSIX" allows me to use the "floor" function below.
                my $Conc = $ARGV[0]; # DNA concentration, ng/µl.
                my $Vol = $ARGV[1];  # DNA volume, µl.
                my $FragMass = $ARGV[2]; # Average DNA fragment mass, ng/fmol.
                my $PropP5 = 0.8058; # Expected proportion of fragment ends that stick to P5 adapters.
                my $PropP7 = 0.1942; # Expected proportion of fragment ends that stick to P7 adapters.
                my $AdapterStockConc = $ARGV[3]; # Annealed adapter stocks molarity (µM, or pmol / µl). Assumed equal for both adapters.
                my $Ends = 2 * $Conc * $Vol / ($FragMass * 1000); # Total number of fragment ends per sample (pmol).
                my $FoldExcess = 10;
                my $Adapter_P5 = $Ends * $PropP5 * $FoldExcess; # Amount of P5 adapter per sample (pmols).
                my $Adapter_P7 = $Ends * $PropP7 * $FoldExcess; # Amount of P7 adapter per sample (pmols).
		my $AdapterWorkVol = 14; # Predetermined volume of adapter working stock, to re-use working stocks across pools (µl).
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

                my $TotVol = ceil(($Vol + 2 + $AdapterVolRxn) / (1 - 1/10 - 1/40));
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
                print "2 µl T4 Ligase, and water to $TotVol µl.\n\n";

                open(LIG, ">>", "ligations2.txt") || die "I cannot open ligations2.txt.\n";
                open(WRK, ">>", "workingStock2.txt") || die "I cannot open workingSotck2.txt.\n";
                my $water = $TotVol - $Vol - $AdapterVolRxn - $TotVol/10 - $TotVol/40 - 2;
                print LIG $ARGV[4], "\t", $ARGV[5], "\t", $ARGV[6], "\t", $Vol, "\t", sprintf("%.2f\t%.2f\t%.2f\t2.0\t%.2f\t%.2f\n", $AdapterVolRxn, $TotVol/10, $TotVol/40, $water, $TotVol);
                print WRK $ARGV[4], "\t", $ARGV[5], "\t", sprintf("%.2f\t%.2f\t%.2f\t%.2f\n", $AdapterStockVol_P5, $AdapterStockVol_P7, $AnnealBufferVol, $AdapterWorkVol);
                close(LIG);
                close(WRK);' ${Concentration[$i]} ${Volume[$i]} ${FragMass[$i]} $StockConc ${Sample[$i]} ${Adapter[$i]} ${Pool[$i]} >> results_2
done
