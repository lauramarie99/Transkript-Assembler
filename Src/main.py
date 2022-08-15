# IMPORT
import sys
import parse_graph_new
import path_enumeration
import pairedbin_enumeration
import networkx as nx
import time
import json
import matplotlib.pyplot as plt
from copy import deepcopy
import optimize
import os

# VARIABLES
start = time.time()
start_gene = time.time()
no_trans = 0
# gene_id = 0
# file_gtf = open("transcripts.gtf", "w")
data_dict = dict()
#read in file to estimate calc time:

# MAIN
if(sys.argv[1] =="-help"):
    print("usage: type python main.py [input graph] [arguments]")
    print("arguments:")
    print("-full for full path enumaration")
    print("-multi for multi bin enumaration")
    print("-paired for paired bin enumaration")
    print("-paired2 for second paired bin enumaration function")
    print("-opt for optimization function and to gain expression levels")
    print("--> requires prior specification of enumaration type")
    print("--> example: main.py Test.graph -paired -opt")
    print("--> results are stored in same folder as save.jsn")
    print("-fullgraph: combine with other arguments to use the full graph (cleaned graph is used otherwise)")
else:
    with open(sys.argv[1]) as file:
        lines = file.readlines()
        num_genes = len([line for line in lines if "--------------" in line])



    with open(sys.argv[1]) as f:
        fileEndReached = False
        f.readline()  # skip ---- seperator line
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
            if ("-fullgraph" not in sys.argv):
                # FULL PATH ENUMERATION
                if (sys.argv[2] =="-full"):

                    transcripts = path_enumeration.enumeration(G_clean,[],"0",["0"],"1",False)
                    transcripts_edge = path_enumeration.enumeration(G_clean, [], "0", ["0"], "1", True)
                    no_trans = no_trans + len(transcripts)

                if (sys.argv[2] == "-multi"):
                # MULTI BIN ENUMERATION
                    multi_bins = path_enumeration.get_multibins(Bins)
                    transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],multi_bins,"1",False)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_clean, [], "0", ["0"], [], multi_bins, "1", True)

                    no_trans = no_trans + len(transcripts)


                # PAIRED BIN ENUMERATION 1
                elif(sys.argv[2] =="-paired"):

                    multi_bins = path_enumeration.get_multibins(Bins)
                    paired_bins = pairedbin_enumeration.get_pairedbins(G_clean,PairedBins_copy,multi_bins)
                    transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],paired_bins+multi_bins,"1",True)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_clean, [], "0", ["0"], [], paired_bins + multi_bins, "1", False)
                    no_trans = no_trans + len(transcripts)


                # PAIRED BIN ENUMERATION 2
                elif(sys.argv[2] == "-paired2"):

                    pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins_copy)
                    multi_bins = path_enumeration.get_multibins(Bins)
                    transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],multi_bins,"1",True)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_clean, [], "0", ["0"], [], multi_bins, "1", False)
                    transcripts_copy = deepcopy(transcripts)
                    filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts_copy,pairedbins_grouped)
                    no_trans = no_trans + len(filtered_transcripts)

                # Optimization
                if("-opt" in sys.argv):
                    var_dict = optimize.model(G_clean, transcripts)
                    data = [(transcripts_edge[i], var_dict[f"expression_levels[{i}]"]) for i in range(len(transcripts_edge))]
                    data_dict[len(data_dict)] = data
                    os.system('clear')
                    print(f"Gene {len(data_dict)} of {num_genes} done.")
                    end_gene = time.time()
                    time_for_gene = end_gene - start_gene
                    start_gene = time.time()
                    print(f"Time for gene: {time_for_gene:.4f}s")
                    time_left = num_genes - len(data_dict) * time_for_gene
                    print(f"Time left:{time_left // 60 // 60:.0f}h {time_left // 60 % 60:.0f}m {time_left % 60:.0f}s")
                    with open("save.json", 'w') as file:
                        json.dump(data_dict, file)

                # ADD TRANSCRIPTS TO GTF FILE
                """
                elif(sys.argv[2] == "-GTF"):
                    gene_id = 1
                    transcript_id = 0
                    for transcript in transcripts:
                        transcript_id += 1
                        parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(gene_id),"Transcript"+str(transcript_id))
                """
            else:
                if (sys.argv[2] =="-full"):

                    transcripts = path_enumeration.enumeration(G_full,[],"0",["0"],"1",False)
                    transcripts_edge = path_enumeration.enumeration(G_full, [], "0", ["0"], "1", True)
                    no_trans = no_trans + len(transcripts)

                if (sys.argv[2] == "-multi"):
                # MULTI BIN ENUMERATION
                    multi_bins = path_enumeration.get_multibins(Bins)
                    transcripts = path_enumeration.enumeration_bins2(G_full,[],"0",["0"],[],multi_bins,"1",False)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_full, [], "0", ["0"], [], multi_bins, "1", True)

                    no_trans = no_trans + len(transcripts)


                # PAIRED BIN ENUMERATION 1
                elif(sys.argv[2] =="-paired"):

                    multi_bins = path_enumeration.get_multibins(Bins)
                    paired_bins = pairedbin_enumeration.get_pairedbins(G_full,PairedBins_copy,multi_bins)
                    transcripts = path_enumeration.enumeration_bins2(G_full,[],"0",["0"],[],paired_bins+multi_bins,"1",True)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_full, [], "0", ["0"], [], paired_bins + multi_bins, "1", False)
                    no_trans = no_trans + len(transcripts)


                # PAIRED BIN ENUMERATION 2
                elif(sys.argv[2] == "-paired2"):

                    pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins_copy)
                    multi_bins = path_enumeration.get_multibins(Bins)
                    transcripts = path_enumeration.enumeration_bins2(G_full,[],"0",["0"],[],multi_bins,"1",True)
                    transcripts_edge = path_enumeration.enumeration_bins2(G_full, [], "0", ["0"], [], multi_bins, "1", False)
                    transcripts_copy = deepcopy(transcripts)
                    filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts_copy,pairedbins_grouped)
                    no_trans = no_trans + len(filtered_transcripts)

                # Optimization
                if("-opt" in sys.argv):
                    var_dict = optimize.model(G_full, transcripts)
                    data = [(transcripts_edge[i], var_dict[f"expression_levels[{i}]"]) for i in range(len(transcripts_edge))]
                    data_dict[len(data_dict)] = data
                    os.system('clear')
                    print(f"Gene {len(data_dict)} of {num_genes} done.")
                    end_gene = time.time()
                    time_for_gene = end_gene - start_gene
                    start_gene = time.time()
                    print(f"Time for gene: {time_for_gene:.4f}s")
                    time_left = num_genes - len(data_dict) * time_for_gene
                    print(f"Time left:{time_left // 60 // 60:.0f}h {time_left // 60 % 60:.0f}m {time_left % 60:.0f}s")
                    with open("save.json", 'w') as file:
                        json.dump(data_dict, file)






        # PRINT RESULTS
        end = time.time()
        print("Gesamtanzahl Transkripte: ", no_trans)
        print('{:5.3f}s'.format(end - start))
        # file_gtf.close()
