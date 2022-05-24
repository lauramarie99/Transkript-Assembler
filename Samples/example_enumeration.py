"""
This is an example for the path enumeration function without any constraints.
The function enumeration() is searching for all possible paths in a given graph.
"""

# IMPORT
from collections import namedtuple
import parse_graph_new
import path_enumeration
import networkx as nx

# MAIN
with open("/home/laura/Documents/Transkript_Assembly/data/human_geuvadis/test3.graph") as f: # Test.graph consists of one single gene entry
    fileEndReached = False
    f.readline()
    while not fileEndReached:
        f.readline() 
        
        # READ META AND BIN DATA
        Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
        Bins = parse_graph_new.parse_bins(f)
        PairedBins = parse_graph_new.parse_pairs(f)
        
        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons)

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()  
            fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons)

        # ENUMERATION
        transcripts = []
        transcripts = path_enumeration.enumeration(G_clean,transcripts,"0",["0"],"1")
        print("Transcripts:", len(transcripts),transcripts)
    
    f.close()
