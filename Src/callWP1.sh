#!/bin/zsh

touch ../Results/resultsWP1.csv
mkdir ../Results/WP1
for fileName in human_geuvadis_simulated_5sets human_leukemia_real_5sets
do
        mkdir ../Results/WP1/$fileName
        graphFileName=$fileName".graph"
        gtfFileName=$fileName".gtf"
        for mode in full multi paired paired2
        do
                mkdir ../Results/WP1/$fileName/$mode
                echo $fileName $mode
                outputFileName=$fileName$mode
                conda activate TranscriptReconstruction
                echo python main.py $graphFileName -$mode -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP1.csv
                python main.py $graphFileName -$mode -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP1.csv
                conda deactivate
                mv $outputFileName".json" ../Results/WP1/$fileName/$mode
                conda activate Cufflinks
                sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP1.csv
                conda deactivate
                mv r.* ../Results/WP1/$fileName/$mode
                mv transcripts.gtf ../Results/WP1/$fileName/$mode/$outputFileName"_transcripts.gtf"
	done
done

gtfFileName=humanDiverse.gtf
for fileName in SRR307903 SRR307911 SRR315323 SRR315334 SRR387661 SRR534291 SRR534307 SRR534319 SRR545695 SRR545723 ERR188021
do
        mkdir ../Results/WP1/$fileName
        graphFileName=$fileName".graph"
        for mode in full multi paired paired2    
        do
                mkdir ../Results/WP1/$fileName/$mode
                echo $fileName $mode   
                outputFileName=$fileName$mode   
                conda activate TranscriptReconstruction
                echo "python main.py $graphFileName -$mode -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP1.csv"
                python main.py $graphFileName -$mode -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP1.csv                                
                conda deactivate
                mv $outputFileName".json" ../Results/WP1/$fileName/$mode
                conda activate Cufflinks
                sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP1.csv
                conda deactivate
                mv r.* ../Results/WP1/$fileName/$mode
                mv transcripts.gtf ../Results/WP1/$fileName/$mode/$outputFileName"_transcripts.gtf"
        done
done
