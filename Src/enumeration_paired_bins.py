# 90841 Transkripte for clean graph

import networkx as nx
import parse_graph_new
import path_enumeration
from collections import namedtuple

BinT = namedtuple('BinT', 'exons count')

# GET_PAIREDBINS FUNCTION
"""
The get_pairedbins function takes all pairedbins and returns a list of bins, which can be used for enumeration.
For each pairedbin, all possible paths connecting the left and right exons are determined.
A new bin is creating for each path, containing the left and right bins and the corresponding path.
"""

def get_pairedbins(graph,pairedbins,multibins):
    all_pairedbins = [] # List storing all bins
    for pairedbin in pairedbins:
        start_node = "" # Start node of enumeration
        end_node = "" # End node of enumeration
        new_pairedbins = [] # New bins, representing all possible paths from the left to the right exons of each paired bin
        left = pairedbin.leftExons # Left exons of paired bin
        right = pairedbin.rightExons # Right exons of paired bin
        count = pairedbin.count # Number of reads found
        
        # CASE 1: All right exons are smaller than the left exons
        if left[0] > right[len(right)-1]:          
            x = left
            left = right
            right = x
        
        # CASE 2: Right and left exons share exons (1,2,4-4,5)
        # CASE 3: Right and left exons do not share exons (1,2-5)
        invalid = False
        # For Loop removes all repeats in the right exons, for example left:1,2,3 and right:3,4,5 will be transformed into left:1,2,3 and right:4,5
        for i in range(0,len(left)):
            if left[i] == right[0]:
                right.pop(0)
                if len(right)!=0:
                    if i!=(len(left)-1) and left[i+1] != right[0]: # A paired bin with left exons 1,2,4 and right exons 2,5,6 is invalid
                        invalid = True
                        break                    
            if len(right) == 0:
                break
        
        if invalid == True:
            continue
        
        
        # If the right exon list is now empty, no enumeration has to be carried out.
        if len(right) == 0:
            new_bin = BinT(exons=left, count=count)
            if ((new_bin not in all_pairedbins) and (len(new_bin.exons) > 2)):
                all_pairedbins.append(new_bin)
            continue
        
        # If the first right exon is smaller than the last left exon, the bin is invalid
        elif right[0] < left[len(left)-1]:
            continue
        
        else:
            start_node = str(left[len(left)-1] + 2) # The start exon is the last exon in the left list. The name of the corresponding node in the graph has to be calculated
            end_node = str(right[0] + 2) # Same for the end exon, which is the first exon in the right exon list
            if (graph.has_node(start_node) and graph.has_node(end_node)):
                new_pairedbins = path_enumeration.enumeration_bins2(graph,new_pairedbins,start_node,[start_node],[],multibins,end_node) # Path enumeration between left and right exons to find all connections
            # Bins have to be completed by adding the missing start and end exons
            for i in range(len(new_pairedbins)):
                if len(left) > 1:
                    exons = left[0:len(left)-1] + new_pairedbins[i] + right
                else:
                    exons = new_pairedbins[i] + right
                # Bin is added to the bin list
                new_bin = BinT(exons=exons, count=count)
                if ((new_bin not in all_pairedbins) and (len(new_bin.exons) > 2)):
                    all_pairedbins.append(new_bin)
    return all_pairedbins
