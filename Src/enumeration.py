# IMPORT
import networkx as nx
from collections import namedtuple
import parse_graph_new

# ENUMERATION FUNCTION
"""
The enumeration function is searching for all possible paths in a given graph.
"""
def enumeration(graph,transcripts:list,node:str,path:list, endnode:str, transcript:bool):
    # if the end node is reached, a path is found and added to the transkript list
    if node == endnode:
        if transcript == True:
            transcripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        else:
            transcripts.append(path)
        return  
    else:
        succ = list(graph.adj[node]) # succ list contains all successors
        for n in succ:
            enumeration(graph,transcripts,n,(path + [n]),endnode,transcript)
    return transcripts


# GET_BINS FUNCTION
"""
The get_bins function returns all bins starting with a specified exon.
The function is needed for the enumeration_bins1 functions. 
"""
def get_bins(bins:list,exon:int):
    new_bins = []
    for bin in bins: 
        if bin.exons[0] == exon:
            new_bins.append(bin)
    return new_bins


# GET_MULTIBINS FUNCTION
"""
The get_multibins function filters all bins with more than two exons. 
The function is needed for the enumeration_bins functions.
"""
def get_multibins(bins:list):
    multi_bins = []
    for bin in bins:
        if len(bin.exons) > 2:
            multi_bins.append(bin)
    return multi_bins


# ENUMERATION_BINS1 FUNCTION
"""
The enumeration_bins1 function is an enumeration function with multibin constraint. 
The bin list contains all bins with more than two exons.
In the beginning act_bins = bins.
"""
def enumeration_bins1(graph,transcripts:list,node:str,path:list,act_bins:list,bins:list,endnode:str,transcript:bool):
    
    if node == endnode:
        if transcript == True:
            transcripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        else:
            transcripts.append(path)
        return
    else:
        succ = list(graph.adj[node]) # succ contains all successor nodes
        for n in succ:
            
            if graph.edges[node,n]['type'] == "SpliceJunction": # check if there's a splice junction: we reach a new exon and have to check for compatibility!
                start_exon = graph.edges[node,n]["startExon"] # start Exon 
                end_exon = graph.edges[node,n]["endExon"] # end Exon
                
                new_bins = [] # stores all active bins
                valid = False # checks if there are compatible bins: Active bins, which contain start and end exon, and the start_exon is not the first exon of the bin.
                other_bins = False # checks if there are active bins, which contain the start exon followed by another exon (not end exon!), and the start_exon is not the first exon of the bin.
                
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            if (i != 0):
                                valid = True
                            if (i != (len(bin1.exons)-2)):
                                new_bins.append(bin1)
                        elif ((start_exon == bin1.exons[i]) and (i != 0)):
                            other_bins = True

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (other_bins == True)):
                    continue

                # Get all new bins starting with the end exon
                new_bins = new_bins + get_bins(bins,end_exon)
            
            # If the edge represents helper edge or exon edge, the bin composition doesn't change
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_bins1(graph,transcripts,n,(path + [n]),new_bins,bins,endnode,transcript)          

    return transcripts

# ENUMERATION_BINS2 FUNCTION
"""
Similar to enumeration_bins1 function, but easier to understand.
The act_bins list is empty at the beginning.
"""
def enumeration_bins2(graph,transcripts:list,node:str,path:list,act_bins:list,bins:list,endnode:str,transcript:bool):
    
    if node == endnode:
        if transcript == True:
            transcripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        else:
            transcripts.append(path)
        return
    else:
        succ = list(graph.adj[node]) # succ contains all successor nodes
        for n in succ:
            
            if graph.edges[node,n]['type'] == "SpliceJunction": # check if there's a splice junction: We reach a new exon and have to check for compatibility!
                start_exon = graph.edges[node,n]["startExon"] # start Exon 
                end_exon = graph.edges[node,n]["endExon"] # end Exon
                
                new_bins = [] # stores all active bins
                valid = False # checks if there are compatible bins: Active bins, which contain start and end exon, and the start_exon is not the first exon of the bin.
                
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            valid = True
                            if (i != (len(bin1.exons)-2)):
                                new_bins.append(bin1)

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (len(act_bins) != 0)):
                    continue

                # Get all new bins
                for bin2 in bins:
                    if ((bin2.exons[0] == start_exon) and (bin2.exons[1] == end_exon)):
                        if len(bin2.exons) > 2:
                            new_bins.append(bin2)
            
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_bins2(graph,transcripts,n,(path + [n]),new_bins,bins,endnode,transcript)          
            
    return transcripts

