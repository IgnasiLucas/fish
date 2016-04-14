Population Genomics of Fish
===========================

I use the directory organization suggested by  William Stafford Noble
(<http://dx.doi.org/10.1371/journal.pcbi.1000424>). On the top level, there
are a few directories with intuitive, logical names: data (not tracked by
git), results, and documents. Below that level, I use a chronological order.
Every day I start an analysis, I name a new folder after the date.

Here, I have a brief summary of each analysis, in reverse chronological order.

2016-04-11
----------
Check if first reads of non-merged pairs cluster together with merged reads.
I am using reads from only one sample, and running vsearch directly. By the
way, vsearch is much more functional than I thougth. Worth learning.


2016-02-23
----------
I run Pyrad's step 3 (clustering within samples) with different values of the
similarity clustering threshold. A high number of Ns requires a lower similarity
clustering threshold. I also plot the distribution of proportion of Ns in reads.
I conclude that the different response of clustering to similarity threshold
between merged and non-merged reads is due to the different proportions of Ns.
Unexpectedly, non-merged, first reads contain higher proportions of Ns.
If many pairs did not merge because of low quality of second reads, as it seems,
it may be worth to try to cluster non-merged, first reads together with merged
reads. A similarity threshold as low as 84% seems adequate. The maximum number
of gaps to include a read in a cluster needs to be raised (default is 3, and it
discards a larger portion of reads, the lower the similarity threshold is).

Note that a more efficient pipeline to determine optimal clustering thresholds
could be implemented using vsearch directly, and/or subsampling the reads.

Conclusion: The clustering threshold should be 86%, for either merged or non-
merged reads. The maximum number of indels, 5.

2016-02-22
----------
I run Pyrad's step 2 (filtering reads) with different settings of the maximum
number of undetermined bases (Ns) to optimize the number of reads recovered.
As expected, merged reads have fewer Ns than non-merged reads. To be able to use
at least 80% of the reads in a sample, I nead to allow at least 65 (non-merged)
or 45 (merged) Ns. Note that this analysis is set to run in the server (64 threads
requested to pyrad).

2016-02-18
----------
For a presentation, I needed some pictures of haplotype networks, and I used
popart. Note that popart is interactive and not called here.

2016-02-16
----------
Here I run SVDQuartets, to infer the species tree from the SNPs. In fact, I am
using the nexus data, including invariable sites and linked SNPs, but apparently
the effect is minimum (Chifman & Kubatko 2014). I save the scores, which tell
me that the power to distinguish among the possible topologies is minimal. However,
the bootstrap agrees that the parallele speciation topology is always the most
supported one.

2016-01-13
----------
Run jmodeltest to determine the best molecular model for each locus, and run
*BEAST with the selected loci.

2016-01-11
----------
Quality check of loci. Select loci for *BEAST analysis and prepare the files.

2015-12-17
----------
Run pyrad twice, with the assembled and the unassembled reads separately.
The run for the unassembled reads failed. Needs to be checked.

2015-12-16b
-----------
I join the two sequencing runs, merge pairs with PEAR, de-multiplex
both assembled and unassembled reads, and finally trim them with
cutadapt. At this point, only unassembled reads may benefit from the
trimming, but PEAR did not accept trimmed reads as input.

2015-12-16
----------
Run pyrad on assembled ('se') and non-assembled ('pe') reads. I note
a large rate of filtered reads, due to low quality. I need to trim the
reads by quality and check for adapters before running pyrad. Deprecated;
see 2015-12-17.

2015-12-15
----------
Use inline barcodes to refine the de-multiplexing done with the indices.
I used SABRE in this step. The reads from the same sample and alternative
runs are joined in the same file. Deprecated; see 2015-12-16b.

2015-12-03
----------
Merging of paired ends using PEAR. A large portion of pairs do not
get assembled. I will run pyrad separately on merged and non-merged
reads. Deprecated; see 2015-12-16b.

2015-12-02
----------
Comparison of the two sequencing runs of the same DNA library. The 
correlation of number of reads per sample is excellent.

