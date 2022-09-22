#!/bin/bash

#call as compareGTF_GTF.sh [Reference.gtf] [Assembly.gtf]
#returns: TP=Number of correctly assembled transcripts  FP=Number of thereby created falsely assembled transcripts

BD=$(dirname "$0")

cuffcompare -r $1 -M -d -e -o r $2 
TP=$(awk '$3 == "=" {print $0}' $BD/r.transcripts.gtf.tmap | wc -l)
P=$(wc -l $BD/r.transcripts.gtf.tmap | cut -d' ' -f4)
FP=$((P-TP))
     
echo -e $d'\t'$TP'\t'$FP

rm r.combined.gtf
rm r.loci
rm r.stats
rm r.tracking
