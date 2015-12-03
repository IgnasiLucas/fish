#/bin/bash
#
#				2015-12-02
#				----------
#
# In order to obtain enough coverage, I sequenced twice the same library, with
# 24 samples in it. Now, it is time to check how the two runs compare. In
# particular, I would like to know:
#
#    1. The effect of the run on sample-specific coverage.
#    2. The effect of the run on sequence quality.
#    3. The effect of the run on sequence error types.
#
# Unfortunately, there is no time to do everything. So, I run only the first
# analysis.
#
# 1. Effect on coverage
#    ------------------
#
DATADIR=`pwd | sed 's/results/data/'`
if [ ! -d $DATADIR ]; then
   mkdir $DATADIR
fi

SAMPLE=(Undetermined St0001 St0003 St0006 St0015 St0016 St0019 St0037 St0039
         St0043 St0044 St0049 St0050 Bl0065 Bl0076 Bl0080 Bl0083
         Bl0104 Bl0108 Bl0091 Bl0093 Bl0094 Bl0095 Bl0098 Bl0116)

for run in 1 2; do
   if [ ! -d $DATADIR/run$run ]; then
      mkdir $DATADIR/run$run
   fi
   for i in `seq 0 24`; do
      if [ ! -e $DATADIR/run$run/${SAMPLE[$i]}_R1.fastq.gz ]; then
         ln -s $DATADIR/../Coregonus$run/${SAMPLE[$i]}_S$i'_L001_R1_001.fastq.gz' $DATADIR/run$run/${SAMPLE[$i]}_R1.fastq.gz
      fi
      if [ ! -e $DATADIR/run$run/${SAMPLE[$i]}_R2.fastq.gz ]; then
         ln -s $DATADIR/../Coregonus$run/${SAMPLE[$i]}_S$i'_L001_R2_001.fastq.gz' $DATADIR/run$run/${SAMPLE[$i]}_R2.fastq.gz
      fi
      if [ ! -e $DATADIR/run$run/${SAMPLE[$i]}_idx.fastq.gz ]; then
         ln -s $DATADIR/../Coregonus$run/${SAMPLE[$i]}_S$i'_L001_I1_001.fastq.gz' $DATADIR/run$run/${SAMPLE[$i]}_idx.fastq.gz
      fi
   done
done

if [ ! -e coverage.txt ] || find -L $DATADIR -type f -newer coverage.txt | grep -q 'fastq'; then
   echo -e "#Sample\tRun\tCoverage" > coverage.txt
   for run in 1 2; do
      for i in `seq 0 24`; do
         COV=` zless $DATADIR/run$run/${SAMPLE[$i]}_idx.fastq.gz | wc -l | gawk '{print $1/4}'`
         echo -e "${SAMPLE[$i]}\t$run\t$COV\t" >> coverage.txt
      done
   done
fi

if [ ! -e coverage.png ] || find . -newer coverage.png | grep -q 'coverage.txt'; then
   if [ ! -e plots01.R ]; then
      echo "plots01.R not found."
   else
      R --save < plots01.R
   fi
fi

# Conclusions
# -----------
#
# After removing the Undetermined reads, the correlation of sample-specific coverage
# between runs is almost perfect. The intercept term is not significant. This is the
# output of the summary of the linear regression of coverage in run 2 on coverage in
# run 1:
#
#   Call:
#   lm(formula = cov[f2] ~ 0 + cov[f1])
#
#   Residuals:
#       Min      1Q  Median      3Q     Max
#   -6590.8 -3148.7  -733.3  1218.3 12955.2
#
#   Coefficients:
#           Estimate Std. Error t value Pr(>|t|)
#   cov[f1]  1.08712    0.00161   675.4   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#   Residual standard error: 4635 on 23 degrees of freedom
#   Multiple R-squared:  0.9999,	Adjusted R-squared:  0.9999
#   F-statistic: 4.562e+05 on 1 and 23 DF,  p-value: < 2.2e-16
#
# This means that, as expected, replicate runs may change the overall mean coverage,
# but they do not modify the distribution of reads among samples, as expected. Note
# that only the change in the number of undetermined reads is not proportional to the
# change of mean coverage (not shown). See the plot "coverage.png" for the correlation.
# Sample St0006 is the one whose sequencing failed, in both runs, perhaps due to
# insufficient starting material or library error.
