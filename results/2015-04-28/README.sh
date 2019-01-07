#!/bin/bash
#
#				2015-04-28
#				----------
#

DATADIR=/home/ilucas/Documents/Coregonus/data/Gagnaire2013

if [ ! -d fastq ]; then
	mkdir fastq
	ln -s $DATADIR/* fastq/
fi

if [ ! -e params.txt ]; then
	pyrad -n
	sed -i '/## 6. /c\TGCA                      ## 6. Restriction overhang.' params.txt
fi

if [ ! -d edits ]; then
	pyrad -p params.txt -s 2
fi
