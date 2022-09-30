# IMPORT
import networkx as nx
from collections import namedtuple
import parse_graph_new

# ENUMERATION FUNCTION
"""
The enumeration function is searching for all possible paths in a given graph.
"""
def enumeration(graph,transcripts:list,node:str,path:list, endnode:str, transcript:bool, maxTranscripts):
    if len(transcripts)>maxTranscripts:
        return
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
            enumeration(graph,transcripts,n,(path + [n]),endnode,transcript, maxTranscripts)
    return transcripts


# GET_MULTIBINS FUNCTION
"""
The get_multibins function returns all bins with more than two exons. 
"""
def get_multibins(bins:list):
    multi_bins = []
    for bin in bins:
        if len(bin.exons) > 2:
            multi_bins.append(bin)
    return multi_bins

  
# ENUMERATION_BINS2 FUNCTION
"""
Enumeration function with multi bin constraint.
The act_bins list is empty at the beginning.
"""
def enumeration_bins2(graph,transcripts:list,node:str,path:list,act_bins:list,bins:list,endnode:str,transcript:bool, maxTranscripts, invalidPathCounter):
    if len(transcripts)> int(1e4) and transcript == True:
        return
    if len(transcripts)>maxTranscripts:
        return
    if invalidPathCounter[0]>1e3:
        return
    if node == endnode:
        if transcript == True:
            transcripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        else:
            transcripts.append(path)
        return
    else:
        succ = list(graph.adj[node]) # succ contains all successor nodes
        for n in succ:
            if n=='1' and endnode!='1':
                invalidPathCounter[0]= invalidPathCounter[0] +1 
                continue
            if graph.edges[node,n]['type'] == "SpliceJunction": # if the edge is a splice junction, we will reach a new exon and have to check for compatibility!
                start_exon = graph.edges[node,n]["startExon"] # start Exon 
                end_exon = graph.edges[node,n]["endExon"] # end Exon
                
                new_bins = [] # stores all active bins
                valid = False # checks if there are compatible bins: Active bins, which contain start and end exon
                
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            valid = True
                            if (i != (len(bin1.exons)-2)):
                                new_bins.append(bin1)

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (len(act_bins) != 0)):
                    invalidPathCounter[0]= invalidPathCounter[0] +1 
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
            try:
                enumeration_bins2(graph,transcripts,n,(path + [n]),new_bins,bins,endnode,transcript, maxTranscripts, invalidPathCounter)          
            except RecursionError as re:
                print('Transkript/Bins exceeds max RecursionDepth')
                return
    return transcripts
