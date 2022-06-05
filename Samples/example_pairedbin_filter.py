"""
This is an example for the paired bins post filter 
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
with open("/home/laura/Documents/Transkript_Assembly/data/human_geuvadis/test7.graph") as f: # Test.graph consists of one single gene entry
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
        multi_bins = path_enumeration.get_multibins(Bins)
        pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins)
        """
        for key in pairedbins_grouped:
            print(key,pairedbins_grouped[key])
        """
        
        transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],multi_bins,"1")
        print("Transkripte vor Filtern:", transcripts)
        filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts,pairedbins_grouped)
        #print(pairedbins_grouped)
        
        print("Transkripte nach Filtern:", filtered_transcripts)
        

f.close()