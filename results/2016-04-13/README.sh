#!/bin/bash
#
#				2016-04-13
#				----------
#
# It is clear that some, if not most, non-merged read pairs failed to
# merge not because they come from long templates, but because of the
# poor quality of the bases near the end of the read.
#
# I want to combine merged and non-merged, first reads in a pyrad run.
# The idea is that some non-merged read pairs failed to merge due to
# low quality sequence tracks that prevented the identification of the
# overlap between first and second reads. Thus, at least the first read
# of such pairs could be clustered with its correctly merged peers.
#
# This way, I would increase the coverage at some loci, and reduce the
# duplicity between merged and non-merged clusters. However, discarding
# second reads from non-merged pairs removes information from the analysis.
# It is important to quantify the loss, and eventually try to recover it.
#
# If the names of non-merged, first reads are marked, I will be able to
# identify clusters composed of only non-merged, first reads. Those
# clusters should be formed with the concatenation of first and second
# reads.
#
# Note: I just noted that PEAR reversed the sense of second reads to
# be in the same strand as first reads. Pyrad assumes that this was the
# case if (and only if) the string '.forward' is included in the name
# of the infile. Since I changed those names, back in 2015-12-16b, pyrad
# has been reverse-complementing my second reads erroneously. That is why
# the concatenated reads in pe/edits/*.edit files show large numbers of
# Ns towards the end of the sequence, instead of showing them near the
# first read. I went back to 2015-12-16b and reversed-complemented the
# non-merged second reads back to their original orientation.
#
# Before trying to use in the same analysis merged and non-merged-concatenated
# reads, which can be a headache, it is worth trying an alternative, more
# sensitive merging. Permisive merging is expected to increase the number of
# merged reads. I will use a more sensitive and faster software, called
# bbmerge, from the bbmap package.

