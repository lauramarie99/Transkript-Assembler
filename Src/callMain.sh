#!/bin/zsh

# Loop from 0 to 5

touch ../Results/results.csv

echo "Data'\t'Graph'\t'full'\t'multi'\t'paired'\t'opt'\t'Norm'\t'Sparsity Constraint'\t'Lambda'\t'Mu'\t'flowOptimization'\t'CostFunctionIndex'\t'maxAdditionalEdgeCount'\t'Mode of Backtrack'\t'Name of gtfFile'\t'Name of csv-Resultfile'\t'Name of jsonFile'\t'Time'\t'Number of Transcripts without Optimization'\t'Number of Transcripts with Optimization'\t'Number of genes with 0 Transcripts with Optimization'\t'Average transcript size without Optimization'\t'Standard deviation of transcript size without Optimization'\t'Average transcript size with Optimization'\t'Standard deviation of transcriptSize with Optimization'\t'Number of single Exon transcripts without Optimization'\t'Number of single Exon transcripts with Optimization'\t'True positives'\t'False positives'\t'Total positives'\t'False negatives'\t'Total Transcripts of ReferenceGTF'\t'Sensitivity on IntronChainLevel'\t'Precision on IntronChainLevel'\t'Fuzzy Sensitivity on IntronChainLevel'\t'Fuzzy Precision on IntronChainLevel" >> ../Results/results.csv
source callWP1.sh
source callWP3.sh
source callWP2.sh
