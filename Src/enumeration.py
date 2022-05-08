
# ENUMERATION FUNCTION
"""
The enumeration function is searching for all possible paths in a given graph
"""
def enumeration(graph,transkripts:list,node:str,path:list):
    # if the end node is reached, a path is found and added to the transkript list
    if node == "1":
        transkripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        return
    else:
        succ = list(graph.adj[node]) # succ list contains all successors
        for n in succ:
            enumeration(graph,transkripts,n,(path + [n]))
    return transkripts
 

# GET BIN FUNCTION
"""
The get_bins function returns all bins starting with a specified exon.
The function is needed for the enumeration_bins functions. 
"""
def get_bins(bins:list,exon:int):
    new_bins = []
    for bin in bins: 
        if bin.exons[0] == exon:
            new_bins.append(bin)
    return new_bins
  
  
# GET MULTI-BINS FUNCTION
"""
The get_multibins function filters all bins with more than one exon. The function is used for the enumeration_bins function
"""
def get_multibins(bins:list):
    multi_bins = []
    for bin in bins:
        if len(bin.exons) > 1: # 2
            multi_bins.append(bin)
    return multi_bins
  
  
# ENUMERATION_BINS FUNCTION
"""
The enumeration_bins function considers only paths with corresponding bins. 
The bins list contains all bins with more than one exon. 
In the beginning all bins (> 1 exon) are handed over (act_bins=bins)
"""
def enumeration_bins(graph,transkripts:list,node:str,path:list,act_bins:list,bins:list):
    if node == "1":
        transkripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        return
    else:
        succ = list(graph.adj[node])
        for n in succ:
            # Change active bin composition if there is a splice junction
            if graph.edges[node,n]['type'] == "SpliceJunction":
                start_exon = graph.edges[node,n]["startExon"]
                end_exon = graph.edges[node,n]["endExon"]
                
                # Check if there are compatible bins
                valid_bins = []
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            valid_bins.append(bin1)
                
                # If there are no compatible bins, stop execution of this path
                if len(valid_bins) == 0:
                    continue
                
                # Get all new bins
                start_bins = get_bins(bins,end_exon)
                
                # Remove all inactive bins
                for bin2 in valid_bins:
                    if bin2.exons[len(bin2.exons)-1] == end_exon:
                        valid_bins.remove(bin2)
                
                # New bins consists of all active valid bins and new bins
                new_bins = start_bins + valid_bins
                
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # follow the path
            enumeration_bins(graph,transkripts,n,(path + [n]),new_bins,bins)          
            
    return transkripts
