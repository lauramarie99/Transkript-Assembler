"""
This is an example for the enumeration function with bin constraint. 
As input file, a graph file with one single gene is used.
The function enumeration_bins considers paths with corresponding multi-bins.
"""

# IMPORT
from collections import namedtuple
from pydoc import pathdirs
import parse_graph_new
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
        print("Bins:", Bins)
        PairedBins = parse_graph_new.parse_pairs(f)
        
        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons) # Full Graph
        

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()  
            fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons) # Cleaned Graph
            nx.draw_networkx(G_full, with_labels=True, arrowsize=12)
            plt.show()
            

        # PATH ENUMERATION OF GRAPH

        transkripts_bins1 = [] # All transkripts obtained by enumeration function 1
        transkripts_bins2 = [] # All transkripts obtained by enumeration function 2

        multi_bins = path_enumeration.get_multibins(Bins) # Filter all bins with more than two exons
        
        transkripts_bins1 = path_enumeration.enumeration_bins1(G_clean,transkripts_bins1,"0",["0"],multi_bins,multi_bins,"1") # Function 1
        print("Transkripts for function 1:",len(transkripts_bins1),transkripts_bins1)

        transkripts_bins2 = path_enumeration.enumeration_bins2(G_clean,transkripts_bins2,"0",["0"],[],multi_bins,"1") # Function 2
        print("Transkripts for function 2:",len(transkripts_bins2),transkripts_bins2)
        
        # The results should be the same!
        

f.close()
