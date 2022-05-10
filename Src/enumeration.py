
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
def get_bins(bins:list,new_bins:list,exon:int):
    for bin in bins: 
        if bin.exons[0] == exon:
            new_bins.append(bin)
    return new_bins
  
  
# GET MULTI-BINS FUNCTION
"""
The get_multibins function filters all bins with more than two exons. The function is used for the enumeration_bins function
"""
def get_multibins(bins:list):
    multi_bins = []
    for bin in bins:
        if len(bin.exons) > 2: # 1
            multi_bins.append(bin)
    return multi_bins
  
  
# ENUMERATION_BINS FUNCTION
"""
The enumeration_bins function considers only paths with corresponding bins. 
The bins list contains all bins with more than two exons. 
In the beginning all bins (> 2 exons) are handed over (act_bins=bins)
"""
def enumeration_bins(graph,transkripts:list,node:str,path:list,act_bins:list,bins:list):
    
    if node == "1":
        transkripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        return
    else:
        succ = list(graph.adj[node]) # succ contains all successor nodes
        for n in succ:
            
            if graph.edges[node,n]['type'] == "SpliceJunction": # check if there's a splice junction (we would reach a new exon and have to check for compatibility!)
                start_exon = graph.edges[node,n]["startExon"] # start Exon 
                end_exon = graph.edges[node,n]["endExon"] # end Exon
                
                new_bins = [] # stores all active bins, which contain the start node followed by the end node and at least one more exon 
                valid = False # checks if there are compatible bins: Active bins, which contain start and end exon, and the start_exon is not the first exon of the bin.
                other_bins = False # checks if there are active bins, which contain the start exon followed by another exon (not end exon!), and the start_exon is not the first exon of the bin.
                
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            if (i != 0):
                                valid = True
                            if (i != len(bin1.exons)-2):
                                new_bins.append(bin1)
                        elif ((start_exon == bin1.exons[i]) and (i != 0)):
                            other_bins = True

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (other_bins == True)):
                    continue

                # Get all new bins starting with the end exon
                new_bins = get_bins(bins,new_bins,end_exon)
            
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_bins(graph,transkripts,n,(path + [n]),new_bins,bins)          
            
    return transkripts
