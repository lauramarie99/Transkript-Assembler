# GET_PAIREDBINS FUNCTION
"""
The get_pairedbins function takes all pairedbins and returns a list of bins, which can be used for enumeration.
For each pairedbin, all possible paths connecting the left and right exons are determined.
A new bin is creating for each path, containing the left and right bins and the corresponding path.
"""
def get_pairedbins(graph,pairedbins):
    all_pairedbins = [] # List storing all bins
    for pairedbin in pairedbins:
        start_node = "" # Start node of enumeration
        end_node = "" # End node of enumeration
        new_pairedbins = [] # New bins, representing all possible paths from the left to the right exons of each paired bin
        left = pairedbin.leftExons # Left exons of paired bin
        right = pairedbin.rightExons # Right exons of paired bin
        
        if left[0] > right[len(right)-1]: # If all right exons are smaller than the left exons, right exons become left exons and vice versa            
            x = left
            left = right
            right = x
        # For Loop removes all repeats in the right exons, for example left:1,2,3 and right:3,4,5 will be transformed into left:1,2,3 and right:4,5
        for i in range(0,len(left)):
            if left[i] == right[0]:
                right.pop(0)                      
            if len(right) == 0:
                break
        # If the right exon list is now empty, no enumeration has to be carried out.
        if ((len(right) == 0) or (left[len(left)-1] > right[0])):
            if ((left not in all_pairedbins) and (len(left) > 2)):
                all_pairedbins.append(left)
        else:
            start_node = str(left[len(left)-1] + 2) # The start exon is the last exon in the left list. The name of the corresponding node in the graph has to be calculated
            end_node = str(right[0] + 2) # Same for the end exon, which is the first exon in the right exon list
            if (graph.has_node(start_node) and graph.has_node(end_node)):
                new_pairedbins = enumeration_between_nodes(graph,new_pairedbins,start_node,[start_node],end_node) # Full Path enumeration between left and right exons to find all connections
            # Bins have to be completed by adding the missing start and end exons
            for i in range(len(new_pairedbins)):
                if len(left) > 1:
                    new_pairedbins[i] = left[0:len(left)-1] + new_pairedbins[i] + right
                else:
                    new_pairedbins[i] = new_pairedbins[i] + right
                # Bin is added to the bin list
                if ((new_pairedbins[i] not in all_pairedbins) and (len(new_pairedbins[i]) > 2)):
                    all_pairedbins.append(new_pairedbins[i])
    return all_pairedbins

# ENUMERATION_BETWEEN_NODES FUNCTION
"""
The enumeration_between_nodes function determines all possible paths between two exons.
"""
def enumeration_between_nodes(graph,transkripts:list,node:str,path:list,endnode):
    # if the end node is reached, a path is found and added to the transkript list
    if node == endnode:
        transkripts.append(parse_graph_new.nodepath_to_transcript(graph,path))
        return
    elif node == "1":
        return
    else:
        succ = list(graph.adj[node]) # succ list contains all successors
        for n in succ:
            enumeration_between_nodes(graph,transkripts,n,(path + [n]),endnode)
    return transkripts

  
# ENUMERATION_PAIREDBINS FUNCTION
"""
The enumeration_pairedbins function is similar to the enumeration_bins2 function.
It determines all possible paths in the given graph using the paired bins.
"""
def enumeration_pairedbins(graph,transkripts:list,node:str,path:list,act_bins:list,bins:list):
    
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
                    for i in range(0,(len(bin1)-1)):
                        if (start_exon == bin1[i]) and (end_exon == bin1[i+1]):
                            valid = True
                            if (i != (len(bin1)-2)):
                                new_bins.append(bin1)

                # If there are no compatible bins, but other active multi-bins, stop execution of this path
                if ((valid == False) and (len(act_bins) != 0)):
                    continue

                # Get all new bins
                for bin2 in bins:
                    if ((bin2[0] == start_exon) and (bin2[1] == end_exon)):
                        if len(bin2) > 2:
                            new_bins.append(bin2)
            
            # If the edge represents helper edge or exon edge, the bin composition is not changed
            else:
                new_bins = act_bins

            # Follow the path
            enumeration_pairedbins(graph,transkripts,n,(path + [n]),new_bins,bins)          
            
    return transkripts


