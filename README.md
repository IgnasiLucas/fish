Population Genomics of Fish
===========================

This is the first project I post on my new Open Science Notebook. I am
still reluctant to share all the details and the data. But I plan to 
open up along the project.

2015-12-02
----------
Comparison of the two sequencing runs of the same DNA library. The 
correlation of number of reads per sample is excellent.

2015-12-03
----------
Merging of paired ends using PEAR. A large portion of pairs do not
get assembled. I will run pyrad separately on merged and non-merged
reads. Deprecated; see 2015-12-16b.

2015-12-15
----------
Use inline barcodes to refine the de-multiplexing done with the indices.
I used SABRE in this step. The reads from the same sample and alternative
runs are joined in the same file. Deprecated; see 2015-12-16b.

2015-12-16
----------
Run pyrad on assembled ('se') and non-assembled ('pe') reads. I note
a large rate of filtered reads, due to low quality. I need to trim the
reads by quality and check for adapters before running pyrad. Deprecated;
see 2015-12-17.

2015-12-16b
-----------
I join the two sequencing runs, merge pairs with PEAR, de-multiplex
both assembled and unassembled reads, and finally trim them with
cutadapt. At this point, only unassembled reads may benefit from the
trimming, but PEAR did not accept trimmed reads as input.

2015-12-17
----------
Run pyrad twice, with the assembled and the unassembled reads separately.
The run for the unassembled reads failed. Needs to be checked.

2016-01-11
----------
Quality check of loci. Select loci for *BEAST analysis and prepare the files.

2016-01-13
----------
Run jmodeltest to determine the best molecular model for each locus, and run
*BEAST with the selected loci.

2016-02-16
----------
Here I run SVDQuartets, to infer the species tree from the SNPs. In fact, I am
using the nexus data, including invariable sites and linked SNPs, but apparently
the effect is minimum (Chifman & Kubatko 2014). I save the scores, which tell
me that the power to distinguish among the possible topologies is minimal. However,
the bootstrap agrees that the parallele speciation topology is always the most
supported one.

2016-02-18
----------
For a presentation, I needed some pictures of haplotype networks, and I used
popart. Note that popart is interactive and not called here.


