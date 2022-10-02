#!/bin/bash

touch ../Results/resultsWP2.csv
datafolder=../Data
for fileName in human_geuvadis_simulated_5sets human_leukemia_real_5sets
do
	mkdir ../Results/WP2/$fileName
    graphFileName=${fileName}.graph
    gtfFileName=${fileName}_ref.gtf
    for mode in full 
    do
		mkdir ../Results/WP2/$fileName/$mode
		for norm in norm1 norm0 norm2
		do
			mkdir ../Results/WP2/$fileName/$mode/$norm
			
			# Write Files for noSparsityConstraint first
			mkdir ../Results/WP2/$fileName/$mode/$norm/NoConstr
			outputFileName=${fileName}_full_${mode}_${norm}_noconstr_transcripts.gtf
			conda activate project1
			echo $fileName $mode $norm 
			echo python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -resultsFilename ../Results/resultsWP2.csv
            python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -resultsFilename ../Results/resultsWP2.csv
            conda deactivate
			conda activate cufflinks
            sh compareGTF_GTF.sh ${datafolder}/*/${gtfFileName} transcripts.gtf >> ../Results/resultsWP2.csv
			conda deactivate
            mv r.* ../Results/WP2/$fileName/$mode/$norm/NoConstr
            mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/NoConstr/${outputFileName}
		
			for constraint in constr1 constr0
			do 
				mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint
				for factor in 0.01 0.1 1 5 10
				do
					mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor
					outputFileName=${fileName}_${mode}_${norm}_${constraint}_${factor}_transcripts.gtf
                    conda activate project1
					echo $fileName $mode $norm $constraint $factor
					echo python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -$constraint -factor $factor-resultsFilename ../Results/resultsWP2.csv
					python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -$constraint -factor $factor -resultsFilename ../Results/resultsWP2.csv
                    conda deactivate
                    conda activate cufflinks
                    sh compareGTF_GTF.sh ${datafolder}/*/${gtfFileName} transcripts.gtf >> ../Results/resultsWP2.csv
                    conda deactivate
                    mv r.* ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor
                    mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor/$outputFileName
                done	
			done
		done	
	done
done


gtfFileName=human_diverse_ref.gtf
for fileName in SRR307903 SRR307911 SRR315323 SRR315334 SRR387661 SRR534291 SRR534307 SRR534319 SRR545695 SRR545723 ERR188021
do
    mkdir ../Results/WP2/$fileName
    graphFileName=${fileName}.graph
    for mode in full
    do
        mkdir ../Results/WP2/$fileName/$mode
        for norm in norm1 norm0 norm2
        do
            mkdir ../Results/WP2/$fileName/$mode/$norm          
            # Write Files for noSparsityConstraint first
            mkdir ../Results/WP2/$fileName/$mode/$norm/NoConstr
            outputFileName=${fileName}_${mode}_${norm}_noconstr_transcripts.gtf
            conda activate project1
            echo $fileName $mode $norm
			echo python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -resultsFilename ../Results/resultsWP2.csv 
			python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -resultsFilename ../Results/resultsWP2.csv
            conda deactivate
            conda activate cufflinks
            sh compareGTF_GTF.sh ${datafolder}/*/${gtfFileName} transcripts.gtf >> ../Results/resultsWP2.csv
            conda deactivate
            mv r.* ../Results/WP2/$fileName/$mode/$norm/NoConstr
            mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/NoConstr/${outputFileName}
                        
            for constraint in constr0 constr1 
            do      
                mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint
                for factor in 0.01 0.1 1 5 10
                do
                    mkdir ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor
                    outputFileName=${fileName}_${mode}_${norm}_${constraint}_${factor}_transcripts.gtf
                    conda activate project1
					echo $fileName $mode $norm $constraint $factor
					echo python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -$constraint -factor $factor -resultsFilename ../Results/resultsWP2.csv
                    python main.py ${datafolder}/*/${graphFileName} -$mode -opt -$norm -$constraint -factor $factor -resultsFilename ../Results/resultsWP2.csv
                    conda deactivate
                    conda activate cufflinks
                    sh compareGTF_GTF.sh ${datafolder}/*/${gtfFileName} transcripts.gtf >> ../Results/resultsWP2.csv
                    conda deactivate
                    mv r.* ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor
                    mv transcripts.gtf ../Results/WP2/$fileName/$mode/$norm/$constraint/$factor/$outputFileName
                done
            done
        done
    done
done