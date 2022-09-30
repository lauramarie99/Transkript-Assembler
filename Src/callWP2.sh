#!/bin/zsh

touch ../Results/resultsWP2.csv
mkdir ../Results/WP2
for fileName in human_geuvadis_simulated_5sets human_leukemia_real_5sets
do
	mkdir ../Results/WP2/$fileName
        graphFileName=$fileName".graph"
        gtfFileName=$fileName".gtf"
        for mode in full
        do
		mkdir ../Results/WP2/$fileName/$mode
		for norm in norm0 norm1
		do
			mkdir ../Results/WP2/$fileName/$mode/$norm
			
			# Write Files for noSparsityConstraint first
			mkdir ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
			outputFileName=$fileName$mode$norm"noSparsityConstraint"
			conda activate TranscriptReconstruction
			echo $fileName $mode $norm 
			echo python main.py $graphFileName -$mode -opt -$norm -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                	python main.py $graphFileName -$mode -opt -$norm -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                	conda deactivate
               	 	mv $outputFileName".json" ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
                	conda activate Cufflinks
                	sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP2.csv
               	 	conda deactivate
                	mv r.* ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
                	mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint/$outputFileName"_transcripts.gtf"	
		
			for constraint in constr0 constr1
			do 
				mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint
				for lambda in 0.01 0.05 0.1 0.5 1 2.5 5 10
				do
					mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
					outputFileName=$fileName$mode$norm$constrain$lambda
                        		conda activate TranscriptReconstruction
					echo $fileName $mode $norm $constraint $lambda
					echo python main.py $graphFileName -$mode -opt -$norm -$constraint $lambda -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
					python main.py $graphFileName -$mode -opt -$norm -$constraint $lambda -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                        		conda deactivate
                        		mv $outputFileName".json" ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
                        		conda activate Cufflinks
                        		sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP2.csv
                        		conda deactivate
                        		mv r.* ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
                        		mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda/$outputFileName"_transcripts.gtf"
                   		done	
			done
		done	
	done
done

gtfFileName=humanDiverse.gtf
for fileName in SRR307903 SRR307911 SRR315323 SRR315334 SRR387661 SRR534291 SRR534307 SRR534319 SRR545695 SRR545723 ERR188021
do
        mkdir ../Results/WP2/$fileName
        graphFileName=$fileName".graph"
        for mode in full
        do
                mkdir ../Results/WP2/$fileName/$mode
                for norm in norm0 norm1
                do
                        mkdir ../Results/WP2/$fileName/$mode/$norm
                        
                        # Write Files for noSparsityConstraint first
                        mkdir ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
                        outputFileName=$fileName$mode$norm"noSparsityConstraint"
                        conda activate TranscriptReconstruction
                        echo $fileName $mode $norm
			echo python main.py $graphFileName -$mode -opt -$norm -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv 
			python main.py $graphFileName -$mode -opt -$norm -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                        conda deactivate
                        mv $outputFileName".json" ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
                        conda activate Cufflinks
                        sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP2.csv
                        conda deactivate
                        mv r.* ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint
                        mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/NoSparsityConstraint/$outputFileName"_transcripts.gtf"
                        
                        for constraint in constr0 constr1 
                        do      
                                mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint
                                for lambda in 0.01 0.05 0.1 0.5 1 2.5 5 10
                                do
                                        mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
                                        outputFileName=$fileName$mode$norm$constrain$lambda
                                        conda activate TranscriptReconstruction
					echo $fileName $mode $norm $constraint $lambda
					echo python main.py $graphFileName -$mode -opt -$norm -$constraint $lambda -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                                        python main.py $graphFileName -$mode -opt -$norm -$constraint $lambda -jsonFilename $outputFileName".json" -resultsFilename ../Results/resultsWP2.csv
                                        conda deactivate
                                        mv $outputFileName".json" ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
                                        conda activate Cufflinks
                                        sh compareGTF_GTF.sh $gtfFileName transcripts.gtf >> ../Results/resultsWP2.csv
                                        conda deactivate
                                        mv r.* ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda
                                        mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/$constraint/$lambda/$outputFileName"_transcripts.gtf"
                                done
                        done
                done
        done
done
