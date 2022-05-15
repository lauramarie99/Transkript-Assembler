
# FULL ENUMERATION FUNCTION
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
The get_multibins function filters all bins with more than two exons. 
The function is needed for the enumeration_bins function
"""
def get_multibins(bins:list):
    multi_bins = []
    for bin in bins:
        if len(bin.exons) > 2:
            multi_bins.append(bin)
    return multi_bins
  
  
# ENUMERATION_BINS FUNCTIONS
"""
The enumeration_bins functions are using the multibin constraint. 
The bins list contains all bins with more than two exons. 
In the beginning act_bins is empty.
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

"""
Similar to enumeration_bins function but easier
"""
def enumeration_bins2(graph,transkripts:list,node:str,path:list,act_bins:list,bins:list):
    
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
                
                for bin1 in act_bins:
                    for i in range(0,(len(bin1.exons)-1)):
                        if (start_exon == bin1.exons[i]) and (end_exon == bin1.exons[i+1]):
                            valid = True
                            if (i != len(bin1.exons)-2):
                                new_bins.append(bin1)

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (len(act_bins) != 0)):
                    continue

                # Get all new bins
                for bin in bins:
                    if bin.exons[0] == start_exon and bin.exons[1] == end_exon:
                        new_bins.append(bin)
            
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_bins2(graph,transkripts,n,(path + [n]),new_bins,bins)          
            
    return transkripts


# GET_PAIREDBINS FUNCTION
"""
The get_pairedbins function takes all paired bins, writes the left and right exons as one bin (exons) and returns an array.
The function is needed for the enumeration_pairedbins function.
"""
def get_pairedbins(pairedbins:list):
    paired_bins = []
    for bin in pairedbins:
        exons = bin.leftExons
        for exon in bin.rightExons:
            if exon not in exons:
                exons.append(exon)
        if ((len(exons) > 1) and (exons not in paired_bins)):
            paired_bins.append(exons)
    return paired_bins

# GET_SUCCEXONS FUNCTION
"""
The get_succesxons function returns all successor exons.
The function is needed for the enumeration_pairedbins function.
"""
def get_succexons(graph,node,succ:list):
    exons = []
    for n in succ:
        if graph.edges[node,n]["type"] == "SpliceJunction":
            exons.append(graph.edges[node,n]["endExon"])
            
    return exons  

# ENUMERATION PAIREDBINS
"""
The enumeration_pairedbins function is using the paired-bin constraint.
The bins list constains all paired bins, which were transformed into single bins by the get_pairedbins function
In the beginning the act_bins list is empty.
"""
def enumeration_pairedbins(graph,transkripts:list,node:str,path:list,act_bins:list,bins:list):
    
    if node == "1":
        transkripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        return
    else:
        succ = list(graph.adj[node]) # succ contains all successor nodes
        succ_exons = get_succexons(graph,node,succ) # get_succexons returns a list containing all successor exons
        
        for n in succ:
            
            if graph.edges[node,n]['type'] == "SpliceJunction": # check if there's a splice junction (we would reach a new exon and have to check for compatibility!)
                start_exon = graph.edges[node,n]["startExon"] # start Exon 
                end_exon = graph.edges[node,n]["endExon"] # end Exon
                
                new_bins = [] # stores all active bins, which contain the end node and at least one more exon 
                valid = False # checks if there are compatible bins: Active bins, which contain the end exon
                other_bins = False # checks if there are active bins which contain another successor (not end_exon)

                for bin1 in act_bins:
                    for i in range(0,len(bin1)):
                        if (end_exon == bin1[i]):
                            valid = True
                            if (i != len(bin1)-1):
                                new_bins.append(bin1)
                        elif (bin1[i] in succ_exons):
                            other_bins = True
                
                # If there are no compatible bins, but active bins containing another successor, stop execution of this path
                if ((valid == False) and (other_bins == True)):
                    continue

                # Get all new bins: Bins containing the start_exon, end_exon and at least one more exon OR bins containing the start_exon, followed by another exon which is not in succ_exons list
                for bin in bins:
                    if (bin[0] == start_exon):
                        if ((bin[1] == end_exon) and (len(bin) > 2)):
                            new_bins.append(bin)
                        elif (bin[1] > end_exon) and (bin[1] not in succ_exons):
                            new_bins.append(bin)
                    
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_pairedbins(graph,transkripts,n,(path + [n]),new_bins,bins)          
            
    return transkripts
