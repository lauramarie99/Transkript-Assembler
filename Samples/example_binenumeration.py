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

        transkripts = [] # All transkripts obtained by "normal" enumeration function
        transkripts_bins = [] # All transkripts obtained by enumeration function with bin constraint

        multi_bins = path_enumeration.get_multibins(Bins) # Filter all bins with more than two exons
        
        transkripts_bins = path_enumeration.enumeration_bins(G_clean,transkripts_bins,"0",["0"],[],multi_bins) # Enumeration with bin-constraint
        print("Transkripts with bin-constraint:",len(transkripts_bins),transkripts_bins)

        transkripts = path_enumeration.enumeration(G_clean,transkripts,"0",["0"]) # Normal enumeration
        print("Transkripts:",len(transkripts),transkripts)
        

f.close()