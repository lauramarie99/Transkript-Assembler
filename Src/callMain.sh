#!/bin/zsh

# Loop from 0 to 5
conda deactivate
graphEnding=.graph
gtfEnding=.gtf
outputFileNameEnding=transcripts.gtf
touch results.txt
pwd=$(pwd)
for fileName in human_geuvadis_simulated_5sets human_leukemia_real_5sets
do
	mkdir $fileName
	graphFileName=$fileName$graphEnding
	gtfFileName=$fileName$gtfEnding
	for value in TLLP TLMF DPLP DPMF
	do
		mkdir $fileName/$value
		for i in $(seq 1 1 5)
		do	
			mkdir $fileName/$value/CostFunction$i
			conda activate TranscriptReconstruction
			python main.py $graphFileName -flowOptimization -$value -costFunction $i -outputFilename transcripts.gtf -additionalEdges 500 -jsonFileName transcripts.json >> results.txt
			conda deactivate
			outputFileName=$fileName$Value$i$outputFileNameEnding
			conda activate Cufflinks
			sh compareGTF_GTF.sh truth.gtf $fileName/$Value/CostFunction$i/transcripts.gtf >> results.txt
			mv transcripts.gtf $fileName/$Value/CostFunction$i/$outputFileName
			mv r.* $fileName/$Value/CostFunction$i/
			conda deactivate
			echo $value$i
		done
	done
done

gtfFileName=humanDiverse.gtf
for graphFileName in SRR307903.graph SRR307911.graph SRR315323.graph SRR315334.graph SRR387661.graph SRR534291.graph SRR534307.graph SRR534319.graph SRR545695.graph SRR545723.graph ERR188021.graph
do
        echo $graphFileName
        for value in TLLP TLMF DPLP DPMF
        do
                for i in $(seq 1 1 5)
                do
                        echo $value$i
                done
        done
done
