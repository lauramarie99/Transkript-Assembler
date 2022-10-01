# Packages
import sys
import parse_graph_new
import path_enumeration
import pairedbin_enumeration
import flowProblem
import networkx as nx
import time
import json
from copy import deepcopy
import optimize
import os
import statistics

# Global variables
start = time.time() # calculate total time needed
start_gene = time.time() # calculate time needed for one gene
no_transcripts = [] # total number of transcripts (without optimization)
no_optimizedTranscripts = [] # list with number of optimized transcripts per gene
numberGenesZeroTranscripts = 0 
numberSingleExonTranscriptsBeforeOptimization = 0 # total number of transcripts with one exon before optimization
numberSingleExonTranscriptsAfterOptimization = 0 # total number of transcripts with one exon after optimization
optimizedTranscriptSize = [] # list with number of exons for each transcript before optimization 
unoptimizedTranscriptSize = [] # list with number of exons for each transcript after optimization
recursionExceededCounter = 0 # number of genes with recursion error
numberGenesFailedOpt = 0 # number of genes with failed optimization
genesFailedOpt = [] # list of genes which failed
data_dict = dict() # dictionary storing all transcripts for each gene
percentageCounter = 0 # percentage of genes completed
transcriptExceededCounter = 0 # number of genes where the number of transcripts exceeds the allowed number (default: 1e4)
maxRecursionDepth = 0 # maximum recursion depth
pairedBinExceeded = 0 # number of genes where the number of paired bins exceeds the allowed number (1e4)
no_Exons = [] # list with number of exons per gene

# Additional Options
# 1. OutputFilename (GTF file)
gtfFilename = "transcripts.gtf" # default name
if "-outputGTF" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-outputGTF":
            gtfFilename = str(sys.argv[i+1])
            break

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

# 4. Result FileName (CSV file)
resultsFilename = "results.csv" # default name
if "-resultsFilename" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-resultsFilename":
            resultsFilename = str(sys.argv[i+1])
            break

file_gtf = open(gtfFilename, "w")

# 4. Json Filename
jsonFilename = "transcripts.json"
if ("-jsonFilename" in sys.argv):
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-jsonFilename":
            jsonFilename = str(sys.argv[i+1])
            break

# # 5. MaximumRecursionDepth
# if "-maxRecursion" in sys.argv:
#     for i in range(len(sys.argv)):
#         if sys.argv[i] == "-maxRecursion":
#             maxRecursionDepth = str(sys.argv[i+1])
#             break
# else:
#     maxRecursionDepth = 1000
# #Set maximum recursion depth
# sys.setrecursionlimit(maxRecursionDepth)

# 5. MaxTranscripts
maxTranscripts = int(1e4) # default value
if "-maxTranscripts" in sys.argv:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-maxTranscripts":
            maxExons = str(sys.argv[i+1])
            break

# 6. Lambda and mu
if ('-opt' not in sys.argv):
    factor=None
    factor=None
else:
    if "-factor" in sys.argv:
        for i in range(len(sys.argv)):
            if sys.argv[i] == "-factor":
                factor = str(sys.argv[i+1])
                break
    else: 
        if("-norm0" in sys.argv and "-constr0" in sys.argv):
            factor = 0.1
        elif ("-norm0" in sys.argv and "-constr1" in sys.argv):
            factor = 0.05
        elif ("-norm1" in sys.argv and "-constr0" in sys.argv):
            factor = 10
        elif ("-norm1" in sys.argv and "-constr1" in sys.argv):
            factor= 5
        elif ("-norm2" in sys.argv and "-constr0" in sys.argv):
            factor = 5
        elif ("-norm2" in sys.argv and "-constr1" in sys.argv):
            factor = 2.5
        elif ("-norm0" in sys.argv):
            factor = None
        elif ("-norm1" in sys.argv):
            factor = None
        elif ("-norm2" in sys.argv):
            factor = None
        else:
            factor = None        


# Main

#prints out usage instructions
if(sys.argv[1] =="-help"):
    print("usage: type python main.py [input graph] [arguments]")
    print("arguments:")
    print("-full for full path enumeration")
    print("-multi for multi bin enumeration")
    print("-paired for paired bin enumeration")
    print("-paired2 for second paired bin enumeration function")
    print("-maxTranscripts: Specify cut-off value, how many Transcripts/Gene (default: 100000)")
    print("-maxRecursion: Specify maximum recursion Depth (Default: 1000)")
    print("-opt for optimization function and to gain expression levels")
    print("--> requires prior specification of enumeration type")
    print("--> specification of norm and sparsity constraint (default values are used otherwise")
    print("-factor: Specify penality size for sparsity constraints")
    print("--> example: main.py Test.graph -paired -opt -norm1 -constr0")
    print("--> results are stored in same folder as save.jsn")
    print("-completegraph: Combine with other arguments to use the full graph (cleaned graph is used otherwise)")
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

# Read in file to estimate calculation time
else:
    with open(sys.argv[1]) as file:
        lines = file.readlines()
        num_genes = len([line for line in lines if "--------------" in line])


# Read in file for parsing
    with open(sys.argv[1]) as f:
        fileEndReached = False
        f.readline()  # skip ---- seperator line
        geneCounter = 0
        residualFlowList = []

        while not fileEndReached:
            # Read meta and bin data from file
            f.readline()  # skip ==META
            Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
            Bins = parse_graph_new.parse_bins(f)
            PairedBins = parse_graph_new.parse_pairs(f)
            PairedBins_copy = deepcopy(PairedBins)
            no_Exons.append(len(Exons))
            # Build graph
            G_full = nx.DiGraph()
            fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons)

            if not fileEndReached and not skip:
                G_clean = nx.DiGraph()
                fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons)
            if skip:
                G_clean = G_full

            # Variables
            transcripts = [] # list of transcripts
            number_optimizedTranscripts = 0 # number of optimized transcripts for -opt argument 
            Graph = None
            if("-completegraph" in sys.argv):
                Graph = G_full
            else:
                Graph = G_clean
            
            # Set maxRecursionDepth Max
            maxRecursionDepth = max(maxRecursionDepth, len(Graph.edges()))
            sys.setrecursionlimit(max(1000,len(Graph.edges())))
            
            # FULL PATH ENUMERATION
            if ("-full" in sys.argv):
                #Skip this gene if it contains too many transcripts
                try:
                    transcripts= path_enumeration.enumeration(Graph,[],"0",["0"],"1",False, maxTranscripts)
                except RecursionError as re:
                    recursionExceededCounter +=1
                    print(str(geneCounter) + ' exceeded maxRecursionNumber.')
                    geneCounter +=1
                    continue

            # MULTI BIN ENUMERATION
            elif ("-multi" in sys.argv):
                multi_bins = path_enumeration.get_multibins(Bins)
                invalidPathCounter = []
                invalidPathCounter.append(0)
                try:
                    transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],multi_bins,"1",False, maxTranscripts, invalidPathCounter)
                except RecursionError as re:
                    recursionExceededCounter +=1
                    print(str(geneCounter) + ' exceeded maxRecursionNumber with.')
                    geneCounter +=1
                    continue

            # PAIRED BIN ENUMERATION 1
            elif("-paired" in sys.argv):
                invalidPathCounter = []
                invalidPathCounter.append(0)
                multi_bins = path_enumeration.get_multibins(Bins)
                paired_bins = pairedbin_enumeration.get_pairedbins(Graph,PairedBins_copy,multi_bins, maxTranscripts)
                if len(paired_bins) > 10000:
                    print('Gene ' + str(geneCounter) + ' exceeded max paired_bins with ' + str(len(paired_bins)))
                    pairedBinExceeded +=1
                    geneCounter +=1
                    continue
                try:
                    transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],paired_bins+multi_bins,"1",False, maxTranscripts, invalidPathCounter)
                except RecursionError as re:
                    recursionExceededCounter +=1
                    print(str(geneCounter) + ' exceeded maxRecursionNumber.')
                    geneCounter +=1
                    continue

            # PAIRED BIN ENUMERATION 2
            elif("-paired2" in sys.argv):
                pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins_copy)
                multi_bins = path_enumeration.get_multibins(Bins)
                invalidPathCounter = []
                invalidPathCounter.append(0)
                try:
                    transcripts = path_enumeration.enumeration_bins2(Graph,[],"0",["0"],[],multi_bins,"1",False, maxTranscripts, invalidPathCounter)
                except RecursionError as re:
                    recursionExceededCounter +=1
                    print(str(geneCounter) + ' exceeded maxRecursionNumber.')
                    geneCounter +=1
                    continue
                transcripts_copy = deepcopy(transcripts)
                filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts_copy,pairedbins_grouped)
            
            # Skip this gene if number of transcripts exceed maxTranscripts
            if len(transcripts)>maxTranscripts:
                transcriptExceededCounter +=1
                print(str(geneCounter) + ' exceeds ' + str(maxTranscripts) + ' transcripts.')
                geneCounter+=1
                continue                

            # Save Number of transcripts for this gene
            no_transcripts.append(len(transcripts))

            # WP2 OPTIMIZATION
            if("-opt" in sys.argv):
                # Get transcripts
                if len(transcripts)==0:
                    transcripts = list(nx.all_simple_paths(Graph, '0', '1')) # Use NetworkX build-in method to obtain all possible paths
                if("-norm0" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L0", sparsity_constr="L0", factor=factor)
                elif ("-norm0" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L0", sparsity_constr="L1", factor=factor)
                elif ("-norm1" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L1", sparsity_constr="L0", factor=factor)
                elif ("-norm1" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L1", sparsity_constr="L1", factor=factor)
                elif ("-norm2" in sys.argv and "-constr0" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L2", sparsity_constr="L0", factor=factor)
                elif ("-norm2" in sys.argv and "-constr1" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L2", sparsity_constr="L1", factor=factor)
                elif ("-norm0" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L0", sparsity_constr=None, factor=0)
                elif ("-norm1" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L1", sparsity_constr=None, factor=0)
                elif ("-norm2" in sys.argv):
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L2", sparsity_constr=None, factor=0)
                else:
                    var_dict = optimize.model(G_clean=Graph, transcripts=transcripts, norm="L1", sparsity_constr=None, factor=0) # if no norm is specified, norm1 is used
                
                if var_dict == None:
                    numberGenesFailedOpt += 1
                    genesFailedOpt.append(geneCounter)
                    geneCounter += 1
                    continue
            

                # function to estimate calculation time
                """
                os.system('clear')
                print(f"Gene {len(data_dict)} of {num_genes} done.")
                end_gene = time.time()
                time_for_gene = end_gene - start_gene
                start_gene = time.time()
                print(f"Time for gene: {time_for_gene:.4f}s")
                time_left = num_genes - len(data_dict) * time_for_gene
                print(f"Time left:{time_left // 60 // 60:.0f}h {time_left // 60 % 60:.0f}m {time_left % 60:.0f}s")
                #save transcripts and their expression levels
                with open("save.json", 'w') as file:
                json.dump(data_dict, file)
                """

            elif "-flowOptimization" in sys.argv:
                
                # Get transcripts
                if ("-TLLP" in sys.argv or "-TLMF" in sys.argv):
                    if len(transcripts)==0:
                        if geneCounter ==8182:
                            print('Trying to establish all paths with NetworkX.')
                # Use NetworkX build-in method to obtain all possible paths
                        transcripts = list(nx.all_simple_paths(Graph, '0', '1'))
                transcriptsCopy = deepcopy(transcripts)

                graphCopy = deepcopy(Graph)
                optimizedTranscripts = []
                        
                #print('CostFunctionIndex = ' + str(costFunctionIndex))
                skipOptimization = False
                
                # Catch infeasible models or models that are unbounded below
                try:
                    g_Star, newGraph, flow = flowProblem.writeGStar(Graph, costFunctionIndex, maxAdditionalEdgeCount, geneCounter)
                    if g_Star == 0:
                        skipOptimization = True    
                except nx.NetworkXUnfeasible or nx.NetworkXUnbounded:
                    skipOptimization = True
                    print('Infeasible Model')
                
                # Execute specified option
                if not skipOptimization:
                    # Get transcripts
                    if "-TLLP" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcriptsCopy, 'longestPath', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                        if residualFlow > 0:
                            residualFlowList.append(geneCounter)
                            print('Residual Flow:' + str(residualFlow) + 'left')
                    elif "-TLMF" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcriptsCopy, 'maximumFlow', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                        if residualFlow > 0:
                            residualFlowList.append(geneCounter)
                            print('Residual Flow:' + str(residualFlow) + 'left')
                    elif "-DPLP" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionDP(newGraph, 'longestPath', flow, geneCounter)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                        if residualFlow > 0:
                            residualFlowList.append(geneCounter)
                            print('Residual Flow:' + str(residualFlow) + 'left')
                    elif "-DPMF" in sys.argv:
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionDP(newGraph, 'maximumFlow', flow, geneCounter)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                        if residualFlow > 0:
                            residualFlowList.append(geneCounter)
                            print('Residual Flow:' + str(residualFlow) + 'left')
                    else: 
                        optimizedGeneTranscripts, residualFlow = flowProblem.flowDecompositionWithTranscriptlist(newGraph, transcripts, 'longestPath', flow)
                        optimizedTranscripts.append(optimizedGeneTranscripts)
                if len(optimizedGeneTranscripts) == 0:
                    print('Gene ' +str(geneCounter) + ' exceeded max length of queue in DP')
                    numberGenesFailedOpt += 1
                    genesFailedOpt.append(geneCounter)
                    geneCounter += 1

            if int(geneCounter/num_genes*100)>= percentageCounter:
                print(f"{percentageCounter}  % finished")
                percentageCounter = percentageCounter+20

            # ADD TRANSCRIPTS TO GTF FILE
            data = []
            if ("-flowOptimization" in sys.argv):
                # Add number to optimized transcripts
                if (len(optimizedGeneTranscripts)) == 0:
                    numberGenesFailedOpt += 1 
                # Save number of optimized transcripts for this Gene if -flowOptimization is used
                no_optimizedTranscripts.append(len(optimizedGeneTranscripts))
                for i in range(len(optimizedGeneTranscripts)):
                    # Write path to transcript
                    transcript = parse_graph_new.nodepath_to_transcript(graphCopy, optimizedGeneTranscripts[i][0])
                    # Save length of optimized transcript
                    optimizedTranscriptSize.append(len(transcript))
                    # Check if transcript is a single Exon transcript
                    if len(transcript) ==1:
                        numberSingleExonTranscriptsAfterOptimization +=1
                        print(geneCounter)
                        print(transcript)
                    # Write GTF-Entry for this transcript
                    parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), "Flow: "+str(optimizedGeneTranscripts[i][1]))
                    #Add transcript and flow as a tuple to data (list)
                    transcriptData = (transcript, optimizedGeneTranscripts[i][1])
                    data.append(transcriptData)
                # Write Data for not-optimized Transcripts
                for i in range(len(transcripts)):
                    transcript = parse_graph_new.nodepath_to_transcript(graphCopy, transcripts[i])
                    unoptimizedTranscriptSize.append(len(transcript))
                    if len(transcript)==1:
                        numberSingleExonTranscriptsBeforeOptimization +=1
            else:
                for i in range(len(transcripts)):
                    transcript = parse_graph_new.nodepath_to_transcript(Graph, transcripts[i])
                    if("-opt" in sys.argv) and var_dict[str(i)] > 0:
                        number_optimizedTranscripts = number_optimizedTranscripts+1
                        optimizedTranscriptSize.append(len(transcript))
                        if len(transcript) ==1:
                            numberSingleExonTranscriptsAfterOptimization +=1
                        parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), 'Expression Level' +str(var_dict[str(i)]))
                        #create list that contains transcripts from all genes and their expression levels. List contains dictionary where key is the gene number (position in file) and values are transcripts and expression level
                        data.append((transcript, var_dict[str(i)]))
                    elif ("-opt" not in sys.argv):
                        parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(geneCounter),str(geneCounter)+"."+str(i+1), 'Expressionlevel: Not determined')
                        data.append(transcript)
                    unoptimizedTranscriptSize.append(len(transcript))
                    if len(transcript)==1:
                        numberSingleExonTranscriptsBeforeOptimization +=1
                # Save number of optimized transcripts for this Gene if -opt is used
                no_optimizedTranscripts.append(number_optimizedTranscripts)
            # Write entry for Gene with list of transcripts and expression levels/flow
            data_dict[geneCounter] = data
            geneCounter = geneCounter + 1
            

# PRINT RESULTS
end = time.time()


# Write Json-File Name
json_object = json.dumps(data_dict)
with open(jsonFilename, 'w') as jsonFile:
    jsonFile.write(json_object)
    jsonFile.close()

# Collect Data for Analysis

# Order
# 0.  InputData (specify name)
# 1.  Graph (cleaned or full)
# 2.  -full (1 = yes, 0=no)
# 3.  -multi (1 = yes, 0=no)
# 4.  -paired (1 = yes, 0no)
# 5   -paired2 (1=yes, 0=np)
# 6.  -maxExons 
# 7.  -maxTranscripts
# 8.  -opt (0 = no, 1 = yes)
# 9.  -norm (0 = L0, 1 = L1, 2 = L2, -1 = none)
# 10.  -constraint (0 = sparsity constraint 0, 1 = sparsity constraint 1, -1 = none)
# 11.  Lambda
# 12.  Mu
# 13.  -flowOptimization (0 = no, 1 = yes)
# 14.  -costFunction (0: f(x) = x, 1: f(x) = x/cov(u,v), 2: f(x) = x/sqrt{cov(u,v)}, 3: f(x) = x^2, 4: f(x) = x^2/cov(u,v), 5: f(x) = x^2*length/cov(u,v), -1: none)
# 15. -maxAdditionalEdges (x = number of additionalEdges for quadraticCostFunction, -1: none)
# 16. Mode of FlowDecomposition (TLLP, TLMF, DPLP, DPMF)
# 17. -outputFilename (outputFilename, default: transcripts.gtf)
# 18. -resultsFilename (resultsFilename, default: results.csv)
# 19. -jsonFilename (name of jsonFilename, -1: none)
# 20. Time
# 21. Number of transcripts determined by pathEnumeration
# 22. Number of transcripts determined by optimization (-opt or flowOptimization)
# 23. Number of genes with 0 Transcripts with Optimization

# 24. Average number of transcripts/gene size without Optimization
# 25. Standard deviation of transcripts/gene without Optimization
# 26. Average number of transcripts/gene size with Optimization
# 27. Standard deviation of transcripts/gene without Optimization

# 28. Average transcript size without Optimization
# 29. Standard deviation of transcriptSize without Optimization
# 30. Average transcript size with Optimization
# 31. Standard deviation of transcriptSize without Optimization

# 32. Number of single Exon transcripts before optimization
# 33. Number of single Exon transcripts after optimization
# 34. Number of total Genes
# 35. Number of skipped Genes
# 36. True positives
# 37. False positives
# 38. Total positives
# 40. False negatives
# 41. Total Transcripts of ReferenceGTF
# 42. Sensitivity on IntronChainLevel
# 43. Precision on IntronChainLevel
# 44. Fuzzy Sensitivity on IntronChainLevel 
# 45. Fuzzy Precision on IntronChainLevel

# Write Dictionary with metaData
metaDataHeader = {}

metaDataHeader[0] = 'Data'
metaDataHeader[1] = 'Graph'
metaDataHeader[2] = 'full'
metaDataHeader[3] = 'multi'
metaDataHeader[4] = 'paired1'
metaDataHeader[5] = 'paired2'
metaDataHeader[6] = 'Number of maximum recursion depth'
metaDataHeader[7] = 'Number of maximum Transcripts/Gene'
metaDataHeader[8] = 'opt'
metaDataHeader[9] = 'Norm'
metaDataHeader[10] = 'Sparsity Constraint'
metaDataHeader[11] = 'Lambda'
metaDataHeader[12] = 'Mu'
metaDataHeader[13] = 'flowOptimization'
metaDataHeader[14] = 'CostFunctionIndex'
metaDataHeader[15] = 'maxAdditionalEdgeCount'
metaDataHeader[16] = 'Mode of Backtrack'
metaDataHeader[17] = 'Name of gtfFile'
metaDataHeader[18] = 'Name of csv-Resultfile'
metaDataHeader[19] = 'Name of jsonFile'
metaDataHeader[20] = 'Time'
metaDataHeader[21] = 'Number of Transcripts without Optimization'
metaDataHeader[22] = 'Number of Transcripts with Optimization'
metaDataHeader[23] = 'Number of genes with 0 Transcripts with Optimization'

metaDataHeader[24] = 'Average number of transcripts/gene without Optimization'
metaDataHeader[25] = 'Standard deviation of transcripts/gene without Optimization'
metaDataHeader[26] = 'Average number of transcripts/gene with Optimization'
metaDataHeader[27] = 'Standard deviation of transcripts/gene with Optimization'

metaDataHeader[28] = 'Average transcript size without Optimization'
metaDataHeader[29] = 'Standard deviation of transcriptSize without Optimization'
metaDataHeader[30] = 'Average transcript size with Optimization'
metaDataHeader[31] = 'Standard deviation of transcriptSize with Optimization'

metaDataHeader[32] = 'Average number of exons/Gene'
metaDataHeader[33] = 'Standard deviation of exons/Gene'

metaDataHeader[34] = 'Number of single Exon transcripts without Optimization'
metaDataHeader[35] = 'Number of single Exon transcripts with Optimization'

metaDataHeader[36] = 'Number of total Genes'
metaDataHeader[37] = 'Number of Genes exceeding maxTranscripts'
metaDataHeader[38] = 'Number of Genes exceeding maxRecursion'
metaDataHeader[39] = 'Number of Genes exceeding maxPairedBins'

metaDataHeader[40] = 'True positives'
metaDataHeader[41] = 'False positives'
metaDataHeader[42] = 'Total positives'
metaDataHeader[43] = 'False negatives'
metaDataHeader[44] = 'Total Transcripts of ReferenceGTF'
metaDataHeader[45] = 'Sensitivity on IntronChainLevel'
metaDataHeader[46] = 'Precision on IntronChainLevel'
metaDataHeader[47] = 'Fuzzy Sensitivity on IntronChainLevel'
metaDataHeader[48] = 'Fuzzy Precision on IntronChainLevel'

metadata = {}

metadata[0] = sys.argv[1]
metadata[1] = 'full' if '-fullgraph' in sys.argv else 'cleaned'
metadata[2] = 1 if '-full' in sys.argv else 0
metadata[3] = 1 if '-multi' in sys.argv else 0
metadata[4] = 1 if '-paired' in sys.argv else 0
metadata[5] = 1 if '-paired2' in sys.argv else 0
metadata[6] = maxRecursionDepth
metadata[7] = maxTranscripts
metadata[8] = 1 if '-opt' in sys.argv else 0
metadata[9] = 0 if '-norm0' in sys.argv else 1 if '-norm1' in sys.argv else 2 if '-norm2' in sys.argv else 1 if '-opt' in sys.argv else -1
metadata[10] = 0 if '-constr0' in sys.argv else 1 if '-constr1' in sys.argv else -1
metadata[11] = factor
metadata[12] = -1
metadata[13] = 1 if '-flowOptimization' in sys.argv else 0 
metadata[14] = costFunctionIndex if '-flowOptimization' in sys.argv else -1
metadata[15] = maxAdditionalEdgeCount if '-flowOptimization' in sys.argv and costFunctionIndex>2 else None
metadata[16] = 'TLLP' if ('-TLLP' in sys.argv and '-flowOptimization' in sys.argv) else 'TLMF' if ('-TLMF' in sys.argv and '-flowOptimization' in sys.argv) else 'DPLP' if ('-DPLP' in sys.argv and '-flowOptimization' in sys.argv) else 'DPMF' if ('-DPMF' in sys.argv and '-flowOptimization' in sys.argv) else 'TLLP' if ('-flowOptimization' in sys.argv) else None
metadata[17] = gtfFilename
metadata[18] = resultsFilename
metadata[19] = jsonFilename if '-jsonFilename' in sys.argv else None
metadata[20] = '{:5.3f}s'.format(end - start)
metadata[21] = sum(no_transcripts)
metadata[22] = sum(no_optimizedTranscripts)
metadata[23] = numberGenesFailedOpt

metadata[24] = statistics.mean(no_transcripts) if len(no_transcripts)>0 else None
metadata[25] = statistics.stdev(no_transcripts) if len(no_transcripts)>0 else None
metadata[26] = statistics.mean(no_optimizedTranscripts) if len(no_optimizedTranscripts)>0 else None
metadata[27] = statistics.stdev(no_optimizedTranscripts) if len(no_optimizedTranscripts)>0 else None

metadata[28] = statistics.mean(unoptimizedTranscriptSize) if len(unoptimizedTranscriptSize)>0 else None
metadata[29] = statistics.stdev(unoptimizedTranscriptSize) if len(unoptimizedTranscriptSize)>0 else None
metadata[30] = statistics.mean(optimizedTranscriptSize) if len(optimizedTranscriptSize)>0 else None
metadata[31] = statistics.stdev(optimizedTranscriptSize) if len(optimizedTranscriptSize)>0 else None

metadata[32] = statistics.mean(no_Exons) if len(no_Exons)>0 else None
metadata[33] = statistics.stdev(no_Exons) if len(no_Exons)>0 else None

metadata[34] = numberSingleExonTranscriptsBeforeOptimization
metadata[35] = numberSingleExonTranscriptsAfterOptimization
metadata[36] = geneCounter
metadata[37] = transcriptExceededCounter
metadata[38] = recursionExceededCounter
metadata[39] = pairedBinExceeded

#Write MetaData to file
metadataFile = open(resultsFilename, 'a')
for i in range(len(metadata)-1):
    metadataFile.write(str(metadata[i]) + "\t")
metadataFile.write(str(metadata[len(metadata)-1]))
metadataFile.close()

print('{:5.3f}s'.format(end - start))


