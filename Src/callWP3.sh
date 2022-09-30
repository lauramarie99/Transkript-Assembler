#!/bin/zsh

touch ../Results/resultsWP3.csv
mkdir ../Results/WP3                                        
for fileName in human_geuvadis_simulated_5sets human_leukemia_real_5sets
do
        mkdir ../Results/WP3/$fileName
        graphFileName=$fileName".graph"
        gtfFileName=$fileName".gtf"
        for value in TLLP TLMF DPLP DPMF
        do
                mkdir ../Results/WP3/$fileName/$value
                for i in $(seq 0 1 5)   
                do
			mkdir ../Results/WP3/$fileName/$value/CostFunction$i
                        if [[ $i -lt 3 ]]
			then
				echo $fileName $value $i
				outputFileName=$fileName$value"CostFunction"$i
				conda activate TranscriptReconstruction
				if [ "$value" = "TLLP" ] || [ "$value" = "TLMF" ]
				then
					echo "python main.py $graphFileName -flowOptimization -full -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv"
                                	python main.py $graphFileName -flowOptimization -full -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv 
                                else
					echo "python main.py $graphFileName -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv"
                                        python main.py $graphFileName -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv
				fi
				conda deactivate
				mv $outputFileName".json" ../Results/WP3/$fileName/$value/CostFunction$i
                                conda activate Cufflinks
                                sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP3.csv
				conda deactivate
                                mv r.* ../Results/WP3/$fileName/$value/CostFunction$i
				mv transcripts.gtf ../Results/WP3/$fileName/$value/CostFunction$i/$outputFileName"_transcripts.gtf"
	
			else
                                for edgeCount in 50 100 150 200 250 
                                do
				        mkdir ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount
					outputFileName=$fileName$value"CostFunction"$i$edgeCount
                        		conda activate TranscriptReconstruction
					if [ "$value" = "TLLP" ] || [ "$value" = "TLMF" ]
					then
						echo "python main.py $graphFileName -full -flowOptimization -$value -costFunction $i transcripts.gtf -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv"
                       	 			python main.py $graphFileName -full -flowOptimization -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv
                        		else
						echo "python main.py $graphFileName -flowOptimization -$value -costFunction $i transcripts.gtf -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv"
                                                python main.py $graphFileName -flowOptimization -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP3.csv
					fi
					conda deactivate
					mv $outputFileName".json" ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount
                        		conda activate Cufflinks
                        		sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP3.csv
					conda deactivate
					mv r* ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount 
					mv transcripts.gtf ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount/$outputFileName"_transcripts.gtf"
				done               
			fi  
                done
        done
done

gtfFileName=humanDiverse.gtf
for fileName in SRR307903 SRR307911 SRR315323 SRR315334 SRR387661 SRR534291 SRR534307 SRR534319 SRR545695 SRR545723 ERR188021
do
        mkdir ../Results/WP3/$fileName
        for value in TLLP TLMF DPLP DPMF
        do
		mkdir ../Results/WP3/$fileName/$value
                for i in $(seq 0 1 5)
                do
                        mkdir ../Results/WP3/$fileName/$value/CostFunction$i
                        if [[ $i -lt 3 ]]
                        then
                                echo $fileName $value $i
                                outputFileName=$fileName$value"CostFunction"$i
                                conda activate TranscriptReconstruction

				if [ "$value" = "TLLP" ] || [ "$value" = "TLMF" ]
                                then
					echo python main.py $fileName".gtf" -full -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
					python main.py $fileName".gtf" -full -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
                                        
                                else
					echo  python main.py $fileName".gtf" -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
					python main.py $fileName".gtf" -flowOptimization -$value -costFunction $i -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
                                fi
                                conda deactivate
                                mv $outputFileName".json" ../Results/WP3/$fileName/$value/CostFunction$i
                                conda activate Cufflinks
                                sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/results.csv
                                conda deactivate
                                mv r.* ../Results/WP3/$fileName/$value/CostFunction$i
                                mv transcripts.gtf ../Results/WP3/$fileName/$value/CostFunction$i/$outputFileName"_transcripts.gtf"

                        else
                                for edgeCount in 50 100 150 200 250
                                do
                                        mkdir ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount
                                        outputFileName=$fileName$value"CostFunction"$i$edgeCount
                                        conda activate TranscriptReconstruction
                                        if [ "$value" = "TLLP" ] || [ "$value" = "TLMF" ]
                                	then
						echo "python main.py $fileName".gtf" -flowOptimization -full -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv"
						python main.py $fileName".gtf" -flowOptimization -full -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
                                        else
						echo "python main.py $fileName".gtf" -flowOptimization -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv"
                                                python main.py $fileName".gtf" -flowOptimization -$value -costFunction $i -additionalEdges $edgeCount -jsonFilename $outputFileName".json" -resultsFilename ../Results/results.csv
					fi
					conda deactivate
                                        mv $outputFileName".json" ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount
                                        conda activate Cufflinks
                                        sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/results.csv
                                        conda deactivate
                                        mv r* ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount
                                        mv transcripts.gtf ../Results/WP3/$fileName/$value/CostFunction$i/$edgeCount/$outputFileName"_transcripts.gtf"
                                done
                        fi
                done
        done
done

