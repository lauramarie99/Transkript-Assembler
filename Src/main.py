# IMPORT
from re import A
import sys
import parse_graph_new
import path_enumeration
import pairedbin_enumeration
import flowProblem
import networkx as nx
import time
import json
import matplotlib.pyplot as plt
from copy import deepcopy
import optimize
import os
import statistics

# VARIABLES
start = time.time()
start_gene = time.time()
no_trans = 0
no_optimizedTranscripts = 0
numberGenesZeroTranscripts = 0
numberSingleExonTranscriptsBeforeOptimization = 0
numberSingleExonTranscriptsAfterOptimization = 0
optimizedTranscriptSize = []
unoptimizedTranscriptSize = []
failed_transcripts = 0
failed_transcripts_ls = []
data_dict = dict()
percentageCounter = 0
# Additional Options

# 1. OutputFilename
if "-outputGTF" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-outputGTF":
            gtfFilename = str(sys.argv[i+1])
            break
else:    
    gtfFilename = "transcripts.gtf"

# 2. CostIndex (determined once)
costFunctionIndex = 1 # default value
if "-costFunction" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-costFunction":
            costFunctionIndex = int(sys.argv[i+1])
            break

# 3. Additional edges
maxAdditionalEdgeCount = 100 # default value
if ("-additionalEdges" in sys.argv):
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-additionalEdges":
            maxAdditionalEdgeCount = int(sys.argv[i+1])

# 4. Result FileName
if "-resultsFilename" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-resultsFilename":
            resultsFilename = str(sys.argv[i+1])
            break
else:    
    resultsFilename = "results.csv"

file_gtf = open(gtfFilename, "w")

# 5. Lambda and mu
if ('-opt' not in sys.argv):
    lambda1=None
    mu = None
else:
    if "-lambda" in sys.argv:
        for i in range(len(sys.argv)):
            if sys.argv[i] == "-lambda":
                lambda1 = str(sys.argv[i+1])
                break
    else: 
        if("-norm0" in sys.argv and "-constr0" in sys.argv):
            lambda1=0.1
        elif ("-norm0" in sys.argv and "-constr1" in sys.argv):
            lambda1 = None
        elif ("-norm1" in sys.argv and "-constr0" in sys.argv):
            lambda1 = 10
        elif ("-norm1" in sys.argv and "-constr1" in sys.argv):
            lambda1=None
        elif ("-norm2" in sys.argv and "-constr0" in sys.argv):
            lambda1 = 5
        elif ("-norm2" in sys.argv and "-constr1" in sys.argv):
            lambda1 = None
        elif ("-norm0" in sys.argv):
            lambda1 = 0
        elif ("-norm1" in sys.argv):
            lambda1 = 0
        elif ("-norm2" in sys.argv):
            lambda1 = None
        else:
            lambda1 = None        

    if "-mu" in sys.argv:
        for i in range(len(sys.argv)):
            if sys.argv[i] == "-mu":
                mu = str(sys.argv[i+1])
                break
    else:
        if("-norm0" in sys.argv and "-constr0" in sys.argv):
            mu = None
        elif ("-norm0" in sys.argv and "-constr1" in sys.argv):
            mu = 0.05
        elif ("-norm1" in sys.argv and "-constr0" in sys.argv):
            mu = None
        elif ("-norm1" in sys.argv and "-constr1" in sys.argv):
            mu = 5
        elif ("-norm2" in sys.argv and "-constr0" in sys.argv):
            mu = None
        elif ("-norm2" in sys.argv and "-constr1" in sys.argv):
            mu = 2.5
        elif ("-norm0" in sys.argv):
            mu = None
        elif ("-norm1" in sys.argv):
            mu = None
        elif ("-norm2" in sys.argv):
            mu = None
        else:
            mu = None

# MAIN
#prints out usage instructions

if(sys.argv[1] =="-help"):
    print("usage: type python main.py [input graph] [arguments]")
    print("arguments:")
    print("-full for full path enumeration")
    print("-multi for multi bin enumeration")
    print("-paired for paired bin enumeration")
    print("-paired2 for second paired bin enumeration function")
    print("-opt for optimization function and to gain expression levels")
    print("--> requires prior specification of enumeration type")
    print("--> specification of norm and sparsity constraint")
    print("-lambda: specify penality size for sparsity constrain 0")
    print("-mu: specify the Maximum number of non-zero transcripts for sparsity contrain 1")
    print("--> example: main.py Test.graph -paired -opt -norm1 -constr0")
    print("--> results are stored in same folder as save.jsn")
    print("-completegraph: combine with other arguments to use the full graph (cleaned graph is used otherwise)")
    print("-flowOptimization: Use flow-based optimization for establishing a set of paths")
    print("--> specify how flow decomposition is being attained:")
    print("-- TLLP: (=TranscriptListLongestPath) in every step the total flow will be reduced by the flow of the longest path of the previously generated set of transcripts")
    print("-- TLMF: (=TranscriptListMaximumFlow) in every step the total flow will be reduced by the flow of the path with the maximumFlow of the previously generated set of transcripts")
    print("--> Caution: For the last two options we recommend using -full, because combining -multi, -paired, -paired2 or -opt with flow-based optimization will probably not accomplish in an optimal solution.")
    print("-- DPLP: (=DynamicProgrammingLongestPath) in every step the total flow will be reduced by the flow of the longest path obtained via dynamic programming (prior establishing the set of paths is not needed for this option)")
    print("-- DPMF: (=DynamicProgrammingMaximumFlow): in every step the total flow will be reduced by the flow of the path with maximumFlow obtained via dynamic programming (prior establishing the set of paths is not needed for this option)")
    print("--> Default Value: TLLP (=Transcript List Longest Path")
    print("-costFunctionX: Specify the costFunction used for flow-based optimization with x (e.g. -flowOptimization 0 for f(x) = x, default value is 1)")
    print("--> 0: f(x) = x" )
    print("--> 1: f(x) = x/cov(u,v)")
    print("--> 2: f(x) = x/sqrt(cov(u,v))")
    print("--> 3: f(x) = x^2 (Approximation of x^2 as ((i+y)*(i+y) - i*i)/y with stepsize y")
    print("--> 4: f(x) = x^2/cov(u,v) (Approximation of x^2 as ((i+y)*(i+y) - i*i)/y with stepsize y")
    print("--> 5: f(x) = x^2*length(u,v)/cov(u,v) (Approximation of x^2 as ((i+y)*(i+y) - i*i)/y with stepsize y")
    print("--> 6: f(x) = x^2 (Exact modelling of x^2 as ((i+y)*(i+y) - i*i)/y with stepsize y=1)")
    print("--> 7: f(x) = x^2 (Approximation of x^2 as x(x-1)/2)")
    print("--> 8: f(x) = x^2/cov(u,v) (Approximation of x^2 as x(x-1)/2)")
    print("-additionalEdges x: Specifiy with x how many additional edges should be added for modelling quadratic cost functions. Be aware that computation time might increases massively.")
    print("-outputGTF: Provide a name for the gtf-File that is written with the transcripts with the ending .gtf (Default: transcripts.gtf")
    print("-resultsFilename: Provide a name where the results are being saved. (Default: results.csv")
    print("-jsonFilename: If you want the transcripts of the genes saved in an additionale json file for further usage, provide a filename in which the transcripts with the expressionlevel or flow will be saved.")
#read in file to estimate calculation time
else:
    with open(sys.argv[1]) as file:
        lines = file.readlines()
        num_genes = len([line for line in lines if "--------------" in line])


#read in file for parsing
    with open(sys.argv[1]) as f:
        fileEndReached = False
        f.readline()  # skip ---- seperator line
        geneCounter = 0
        residualFlowList = []

        while not fileEndReached:
            # READ META AND BIN DATA FROM FILE
            f.readline()  # skip ==META
            Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
            Bins = parse_graph_new.parse_bins(f)
            PairedBins = parse_graph_new.parse_pairs(f)
            PairedBins_copy = deepcopy(PairedBins)

            # BUILD GRAPH
            G_full = nx.DiGraph()
            fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons)

            if not fileEndReached and not skip:
                G_clean = nx.DiGraph()
                fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons)
                # nx.draw_networkx(G_clean, with_labels=True, arrowsize=12)
                # plt.show()
            if skip:
                G_clean = G_full

            transcripts = []
            Graph = None
            
            if("-completegraph" in sys.argv):
                Graph = G_full
            else:
                Graph = G_clean

            # FULL PATH ENUMERATION
            if ("-full" in sys.argv):

                transcripts = path_enumeration.enumeration(Graph,[],"0",["0"],"1",False)
                no_trans = no_trans + len(transcripts)

            # MULTI BIN ENUMERATION
            if ("-multi" in sys.argv):

                multi_bins = path_enumeration.get_multibins(Bins)
                transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],multi_bins,"1",False)
                no_trans = no_trans + len(transcripts)

            # PAIRED BIN ENUMERATION 1
            elif("-paired" in sys.argv):

                multi_bins = path_enumeration.get_multibins(Bins)
                paired_bins = pairedbin_enumeration.get_pairedbins(Graph,PairedBins_copy,multi_bins)
                transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],paired_bins+multi_bins,"1",False)
                no_trans = no_trans + len(transcripts)


            # PAIRED BIN ENUMERATION 2
            elif("-paired2" in sys.argv):

                pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins_copy)
                multi_bins = path_enumeration.get_multibins(Bins)
                transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],multi_bins,"1",False)
                transcripts_copy = deepcopy(transcripts)
                filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts_copy,pairedbins_grouped)
                no_trans = no_trans + len(filtered_transcripts)

            # OPTIMIZATION
            if("-opt" in sys.argv):
                if("-norm0" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L0", sparsity_constr="L0", factor=lambda1)
                elif ("-norm0" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L0", "L1", mu)
                elif ("-norm1" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L1", "L0", lambda1)
                elif ("-norm1" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L1", "L1", mu)
                elif ("-norm2" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L2", "L0", lambda1)
                elif ("-norm2" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L2", "L1", mu)
                elif ("-norm0" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L0", None, lambda1)
                elif ("-norm1" in sys.argv):
                    var_dict = optimize.model(Graph, transcripts, "L1", None, lambda1)
                elif ("-norm2" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L2", sparsity_constr=None, factor=0)
                else:
                    var_dict = optimize.model(Graph, transcripts, "L1", None, 0) # if no norm is specified, norm1 is used
                
                if var_dict == None:
                    failed_transcripts += 1
                    failed_transcripts_ls.append(geneCounter)
                    geneCounter += 1
                    continue
            

                #function to estimate calculation time
                os.system('clear')
                print(f"Gene {len(data_dict)} of {num_genes} done.")
                end_gene = time.time()
                time_for_gene = end_gene - start_gene
                start_gene = time.time()
                print(f"Time for gene: {time_for_gene:.4f}s")
                time_left = num_genes - len(data_dict) * time_for_gene
                print(f"Time left:{time_left // 60 // 60:.0f}h {time_left // 60 % 60:.0f}m {time_left % 60:.0f}s")
                
                """
                #save transcripts and their expression levels
                with open("save.json", 'w') as file:
                json.dump(data_dict, file)
                """

            elif "-flowOptimization" in sys.argv:
                
                graphCopy = deepcopy(Graph)
                optimizedTranscripts = []
                        
                #print('CostFunctionIndex = ' + str(costFunctionIndex))
                skipOptimization = False
                
                # Catch infeasible models or models that are unbounded below
                try:
                    g_Star, newGraph, flow = flowProblem.writeGStar(Graph, costFunctionIndex, maxAdditionalEdgeCount)
                    if g_Star == 0:
                        skipOptimization = True    
                except nx.NetworkXUnfeasible or nx.NetworkXUnbounded:
                    skipOptimization = True
                    print('Infeasible Model')
                
                # Execute specified option
                if not skipOptimization:
                    if "-TLLP" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcripts, 'longestPath', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                    elif "-TLMF" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcripts, 'maximumFlow', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                    elif "-DPLP" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionDP(newGraph, 'longestPath', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                    elif "-DPMF" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionDP(newGraph, 'maximumFlow', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                        if residualFlow > 0:
                            residualFlowList.append(geneCounter)
                            print('Residual Flow:' + str(residualFlow) + 'left')
                    else: 
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcripts, 'longestPath', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                    
                # if int(geneCounter/num_genes*100)>= percentageCounter:
                #     print(f"{percentageCounter}  % finished")
                #     percentageCounter = percentageCounter+ 1 


            # ADD TRANSCRIPTS TO GTF FILE
            data = []
            if ("-flowOptimization" in sys.argv):
                if (len(optimizedGeneTranscripts)) == 0:
                    failed_transcripts += 1 
                no_optimizedTranscripts = no_optimizedTranscripts+len(optimizedGeneTranscripts)
                for i in range(len(optimizedGeneTranscripts)):
                    transcript = parse_graph_new.nodepath_to_transcript(graphCopy, optimizedGeneTranscripts[i][0])
                    optimizedTranscriptSize.append(len(transcript))
                    if len(transcript) ==1:
                        numberSingleExonTranscriptsAfterOptimization +=1
                    parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), "Flow: "+str(optimizedGeneTranscripts[i][1]))
                    #create list that contains transcripts from all genes and their expression levels. List contains dictionary where key is the gene number (position in file) and values are transcripts and expression level
                    transcriptData = (transcript, optimizedGeneTranscripts[i][1])
                    data.append(transcriptData)
            else:
                for i in range(len(transcripts)):
                    transcript = parse_graph_new.nodepath_to_transcript(Graph, transcripts[i])
                    if("-opt" in sys.argv) and var_dict[str(i)] > 0:
                        no_optimizedTranscripts = no_optimizedTranscripts+1
                        optimizedTranscriptSize.append(len(transcript))
                        if len(transcript) ==1:
                            numberSingleExonTranscriptsAfterOptimization +=1
                        parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), 'Expression Level' +str(var_dict[str(i)]))
                        #create list that contains transcripts from all genes and their expression levels. List contains dictionary where key is the gene number (position in file) and values are transcripts and expression level
                        data.append((transcript, var_dict[str(i)]))
                    elif ("-opt" not in sys.argv):
                        unoptimizedTranscriptSize.append(len(transcript))
                        if len(transcript)==1:
                            numberSingleExonTranscriptsBeforeOptimization +=1
                        parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), 'Expressionlevel: Not determined')
                        data.append(transcript)
                #print(transcript)
            data_dict[geneCounter] = data
            geneCounter = geneCounter + 1


# PRINT RESULTS
end = time.time()
# print("Number of transcripts: ", no_trans)

#print('Time: ' + '{:5.3f}s'.format(end - start))

# file_gtf.close()
# print("Optimization failed for ", failed_transcripts, " gene")
# print(failed_transcripts_ls)
# # print(data_dict)
# # print(var_dict)
#print(residualFlowList)

# 4. Json-File Name
if ("-jsonFilename" in sys.argv):
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-jsonFilename":
            jsonFilename = str(sys.argv[i+1])
    json_object = json.dumps(data_dict)
    with open(jsonFilename, 'w') as jsonFile:
        jsonFile.write(json_object)
        jsonFile.write('\n')
        jsonFile.write('Time: ' + str(end - start))
        jsonFile.close()

# Collect Data for Analysis

# Order
# 0.  InputData (specify name)
# 1.  Graph (cleaned or full)
# 2.  -full (1 = yes, 0=None)
# 3.  -multi (1 = yes, 0=None)
# 4.  -paired (1 = yes, 0 = None)
# 5   -paired2 (1=yes, else None)
# 6.  -opt (0 = no, 1 = yes)
# 7.  -norm (0 = L0, 1 = L1, 2 = L2, -1 = none)
# 8.  -constraint (0 = sparsity constraint 0, 1 = sparsity constraint 1, -1 = none)
# 9.  Lambda
# 10.  Mu
# 11.  -flowOptimization (0 = no, 1 = yes)
# 12.  -costFunction (0: f(x) = x, 1: f(x) = x/cov(u,v), 2: f(x) = x/sqrt{cov(u,v)}, 3: f(x) = x^2, 4: f(x) = x^2/cov(u,v), 5: f(x) = x^2*length/cov(u,v), -1: none)
# 13. -maxAdditionalEdges (x = number of additionalEdges for quadraticCostFunction, -1: none)
# 14. Mode of FlowDecomposition (TLLP, TLMF, DPLP, DPMF)
# 15. -outputFilename (outputFilename, default: transcripts.gtf)
# 16. -resultsFilename (resultsFilename, default: results.csv)
# 17. -jsonFilename (name of jsonFilename, -1: none)
# 18. Time
# 19. Number of transcripts determined by pathEnumeration
# 20. Number of transcripts determined by optimization (-opt or flowOptimization)
# 21. Number of genes with 0 Transcripts with Optimization
# 22. Average transcript size without Optimization
# 23. Standard deviation of transcriptSize without Optimization
# 24. Average transcript size with Optimization
# 25. Standard deviation of transcriptSize without Optimization
# 26. Number of single Exon transcripts before optimization
# 27. Number of single Exon transcripts after optimization
# 28. True positives
# 29. False positives
# 30. Total positives
# 31. False negatives
# 32. Total Transcripts of ReferenceGTF
# 33. Sensitivity on IntronChainLevel
# 34. Precision on IntronChainLevel
# 35. Fuzzy Sensitivity on IntronChainLevel 
# 36. Fuzzy Precision on IntronChainLevel

# Write Dictionary with metaData
metaDataHeader = {}

metaDataHeader[0] = 'Data'
metaDataHeader[1] = 'Graph'
metaDataHeader[2] = 'full'
metaDataHeader[3] = 'multi'
metaDataHeader[4] = 'paired1'
metaDataHeader[5] = 'paired1'
metaDataHeader[6] = 'opt'
metaDataHeader[7] = 'Norm'
metaDataHeader[8] = 'Sparsity Constraint'
metaDataHeader[9] = 'Lambda'
metaDataHeader[10] = 'Mu'
metaDataHeader[11] = 'flowOptimization'
metaDataHeader[12] = 'CostFunctionIndex'
metaDataHeader[13] = 'maxAdditionalEdgeCount'
metaDataHeader[14] = 'Mode of Backtrack'
metaDataHeader[15] = 'Name of gtfFile'
metaDataHeader[16] = 'Name of csv-Resultfile'
metaDataHeader[17] = 'Name of jsonFile'
metaDataHeader[18] = 'Time'
metaDataHeader[19] = 'Number of Transcripts without Optimization'
metaDataHeader[20] = 'Number of Transcripts with Optimization'
metaDataHeader[21] = 'Number of genes with 0 Transcripts with Optimization'
metaDataHeader[22] = 'Average transcript size without Optimization'
metaDataHeader[23] = 'Standard deviation of transcriptSize without Optimization'
metaDataHeader[24] = 'Average transcript size with Optimization'
metaDataHeader[25] = 'Standard deviation of transcriptSize with Optimization'
metaDataHeader[26] = 'Number of single Exon transcripts without Optimization'
metaDataHeader[27] = 'Number of single Exon transcripts with Optimization'
metaDataHeader[28] = 'True positives'
metaDataHeader[29] = 'False positives'
metaDataHeader[30] = 'Total positives'
metaDataHeader[31] = 'False negatives'
metaDataHeader[32] = 'Total Transcripts of ReferenceGTF'
metaDataHeader[33] = 'Sensitivity on IntronChainLevel'
metaDataHeader[34] = 'Precision on IntronChainLevel'
metaDataHeader[35] = 'Fuzzy Sensitivity on IntronChainLevel'
metaDataHeader[36] = 'Fuzzy Precision on IntronChainLevel'

metadata = {}

metadata[0] = sys.argv[1]
metadata[1] = 'full' if '-fullgraph' in sys.argv else 'cleaned'
metadata[2] = 1 if '-full' in sys.argv else None
metadata[3] = 1 if '-multi' in sys.argv else None
metadata[4] = 1 if '-paired' in sys.argv else None
metadata[5] = 1 if '-paired2' in sys.argv else None
metadata[6] = 1 if '-opt' in sys.argv else 0
metadata[7] = 0 if '-norm0' in sys.argv else 1 if '-norm1' in sys.argv else 2 if '-norm2' in sys.argv else 1 if '-opt' in sys.argv else None
metadata[8] = 0 if '-constr0' in sys.argv else 1 if '-constr1' in sys.argv else None
metadata[9] = lambda1
metadata[10] = mu
metadata[11] = 1 if '-flowOptimization' in sys.argv else 0 
metadata[12] = costFunctionIndex if '-flowOptimization' in sys.argv else None
metadata[13] = maxAdditionalEdgeCount if '-flowOptimization' in sys.argv and costFunctionIndex>2 else None
metadata[14] = 'TLLP' if ('-TLLP' in sys.argv and '-flowOptimization' in sys.argv) else 'TLMF' if ('-TLMF' in sys.argv and '-flowOptimization' in sys.argv) else 'DPLP' if ('-DPLP' in sys.argv and '-flowOptimization' in sys.argv) else 'DPMF' if ('-DPMF' in sys.argv and '-flowOptimization' in sys.argv) else 'TLLP' if ('-flowOptimization' in sys.argv) else None
metadata[15] = gtfFilename
metadata[16] = resultsFilename
metadata[17] = jsonFilename if 'jsonFilename' in sys.argv else None
metadata[18] = '{:5.3f}s'.format(end - start)
metadata[19] = no_trans
metadata[20] = no_optimizedTranscripts
metadata[21] = failed_transcripts
metadata[22] = statistics.mean(unoptimizedTranscriptSize) if len(unoptimizedTranscriptSize)>0 else None
metadata[23] = statistics.stdev(unoptimizedTranscriptSize) if len(unoptimizedTranscriptSize)>0 else None
metadata[24] = statistics.mean(optimizedTranscriptSize) if len(optimizedTranscriptSize)>0 else None
metadata[25] = statistics.stdev(optimizedTranscriptSize) if len(optimizedTranscriptSize)>0 else None
metadata[26] = numberSingleExonTranscriptsBeforeOptimization
metadata[27] = numberSingleExonTranscriptsAfterOptimization

#Write MetaData to file
metadataFile = open(resultsFilename, 'a')
for i in range(len(metadata)-1):
    metadataFile.write(str(metadata[i]) + "\t")
metadataFile.write(str(metadata[len(metadata)-1]))
metadataFile.close()



