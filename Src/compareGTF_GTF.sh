##!/bin/zsh

#call as compareGTF_GTF.sh [Reference.gtf] [Assembly.gtf]
#returns: TP=Number of correctly assembled transcripts  FP=Number of thereby created falsely assembled transcripts

BD=$(dirname "$0")

cuffcompare -r $1 -M -d -e -o r $2 
TP=$(awk '$3 == "=" {print $0}' r.transcripts.gtf.tmap | wc -l | sed s/' '/''/g)
P=$(wc -l r.transcripts.gtf.tmap | cut -d'r' -f1 | sed s/' '/''/g)
FP=$((P-TP))
TotalReference=$(awk '$3 == "transcript" {print}' $1 | wc -l | sed s/' '/''/g)
FN=$((TotalReference-TP))
Sensitivity=$(head -14 r.stats | tail -1 | cut -d":" -f2 | cut -f2 | sed s/' '/''/g)
Precision=$(head -14 r.stats | tail -1 | cut -d":" -f2 | cut -f3 | sed s/' '/''/g)
FuzzySensitivity=$(head -14 r.stats | tail -1 | cut -d":" -f2 | cut -f4 | sed s/' '/''/g)
FuzzyPrecision=$(head -14 r.stats | tail -1 | cut -d":" -f2 | cut -f5 | sed s/' '/''/g)

echo '\t'$TP'\t'$FP'\t'$P'\t'$FN'\t'$TotalReference'\t'$Sensitivity'\t'$Precision'\t'$FuzzySensitivity'\t'$FuzzyPrecision

rm r.combined.gtf
