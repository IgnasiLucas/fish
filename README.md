Population Genomics of Fish
===========================

I use the directory organization suggested by  William Stafford Noble
(<http://dx.doi.org/10.1371/journal.pcbi.1000424>). On the top level, there
are a few directories with intuitive, logical names: data (not tracked),
results, and documents. Below that level, I use a chronological order.
Every day I start an analysis, I name a new folder after the date.

Here, I have a brief summary of each analysis, in reverse chronological order.


2016-09-19
----------
Here I check the effect of the maximum number of Ns allowed in a read on the
number of clusters produced. The relationship is direct, almost lineal. But it
seems to be caused by a corresponding increase in the number of reads available.

2016-09-14
----------
Here, I repeat the mapping of reads from only one sample with different settings
in order to report all the alternative mappings. I confirm that a large portion
of reads map ambituously due to similarities among clusters. A manual alignment
of consensus sequences of similar clusters makes me think that they are actually
from the same locus but trimmed differently. Actually, I had not noticed before
that pyrad trims way too many reads with the strict setting. 

2016-06-29
----------
In rough agreement with 2016-04-11, around 20% of non-merged reads map to the
consensus sequences of clusters made with merged reads. That is to say that
we cannot expect to improve the coverage of merged clusters by more than 20%
when we map non-merged reads to them. The analysis also reveals a potential
problem with the clusters built with merged reads in 2016-05-10: about 50%
of merged reads map ambiguously to more than one consensus sequence. This
implies that the mapping with bowtie2 is more sensitive than the clustering,
and makes me suspect that the clustering is spuriously splitting some loci,
maybe due to an excess of Ns in the reads.

2016-05-10
----------
Before the process finished in folder 2016-05-03 to build consensus sequences
for merged reads, I start the same process again, with 32 threads, to check
if that is going to be much faster. Indeed, after a couple of days, it was
clear that 32 threads speed up the clustering step by a factor of at least 3.

2016-05-03
----------
I implement here the idea of creating a reference genome with the consensus
of the merged reads from all samples, in order to map later both merged and
non-merged reads against this reference. Only unmapped reads would then be
processed separately with pyrad. The only problem here is that to cluster
merged reads (relatively long) with pyrad (which in turn calls vsearch) is
very slow. The estimated time is more than 3 months. Unfortunately, the process
started running with only 6 threads, due to an error.

2016-04-13
----------
I used bbmerge, from the bbmap package, to attempt more permissive merging
of at least one sample, to check if that removes the problem of split clusters
between merged and non-merged reads. Conclusion: de novo assemblies, based
on k-mers, are not more sensitive than PEAR. This option is discarded.


2016-04-11
----------
Check if first reads of non-merged pairs cluster together with merged reads.
I am using reads from only one sample, and running vsearch directly. By the
way, vsearch is much more functional than I thougth. Worth learning. 

Conclusion: Among the 951073 non-merged (unique) reads from BlCl0091, 219029 (23%)
have their most similar read among merged reads. This is an important
number of reads that can increase the coverage of several loci. However,
the number of non-merged-specific loci is even higher, and should not
be discarded.


2016-02-23
----------
I run Pyrad's step 3 (clustering within samples) with different values of the
similarity clustering threshold. A high number of Ns requires a lower similarity
clustering threshold. I also plot the distribution of proportion of Ns in reads.
I conclude that the different response of clustering to similarity threshold
between merged and non-merged reads is due to the different proportions of Ns.
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

