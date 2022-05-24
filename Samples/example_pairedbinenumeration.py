"""
This is an example for the enumeration function with paired-bin constraint. 
As input file, a graph file with one single gene is used.
The function enumeration_pairedbins considers paths with corresponding paired bins.
"""

# IMPORT
from collections import namedtuple
from pydoc import pathdirs
import parse_graph_new
import pairedbin_enumeration
import path_enumeration
import networkx as nx
from matplotlib import pyplot as plt

# MAIN
with open("/home/laura/Documents/Transkript_Assembly/data/human_geuvadis/test3.graph") as f: # Test.graph consists of one single gene entry
    fileEndReached = False
    f.readline()
    while not fileEndReached:
        
        f.readline()
        Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
        Bins = parse_graph_new.parse_bins(f)
        
        PairedBins = parse_graph_new.parse_pairs(f)
        
        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons) # Full Graph
        

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()  
            fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons) # Cleaned Graph
            nx.draw_networkx(G_clean, with_labels=True, arrowsize=12)
            plt.show()
            

        # PATH ENUMERATION OF GRAPH

        filtered_bins = path_enumeration.get_multibins(Bins)
        pairedbins = pairedbin_enumeration.get_pairedbins(G_clean,PairedBins,filtered_bins)
        transkripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],pairedbins,"1")
        print("Transkripts:", len(transkripts), transkripts)
        

f.close()
