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
        invalid = False # Boolean to check if bin is invalid
        
        # CASE 1: The first right exon is smaller than the first left exon
        if left[0] > right[0]:
            x = left
            left = right
            right = x
        # CASE 2: Right and left exons share exons (1,2,4-4,5)
        # CASE 3: Right and left exons do not share exons (1,2-5)
        # For Loop removes all repeats in the right exons, for example left:1,2,3 and right:3,4,5 will be transformed into left:1,2,3 and right:4,5
        for i in range(0,len(left)):
            if left[i] == right[0]:
                right.pop(0)
                if len(right) > 0:
                    if i!=(len(left)-1) and left[i+1] != right[0]: # A paired bin with left exons 1,2,4 and right exons 2,5,6 is invalid
                        invalid = True
                        break  
            if len(right) == 0:
                break
        
        # If current bin is invalid, continue with the next paired bin
        if invalid == True:
            continue

        # If the right exon list is now empty, no enumeration has to be carried out
        if len(right) == 0:
            new_bin = BinT(exons=left, count=count)
            if (len(new_bin.exons) > 2):
                all_pairedbins.append(new_bin)
            continue
        
        # If the first right exon is smaller than the last left exon, the bin is also invalid
        if right[0] < left[len(left)-1]:
            continue
        
        # The start node and end node for the enumeration has to be determined
        start_exon = left[len(left)-1]
        end_exon = right[0]
        for edgeKey, edgeValue in graph.edges.items():
            if edgeValue['type'] == 'Exon':
                if edgeValue['exon'] == start_exon:
                    start_node = edgeKey[0]
                elif edgeValue['exon'] == end_exon:
                    end_node = edgeKey[0]
            if start_node != "" and end_node != "":
                break

        # If start exon or end exon is not found, continue with next paired bin     
        if start_node == "" or end_node == "":
            continue

        # Path enumeration between left and right exons is carried out to find all possible connections
        new_pairedbins = path_enumeration.enumeration_bins2(graph,new_pairedbins,start_node,[start_node],[],multibins,end_node)
            
        # Bins have to be completed by adding the missing start and end exons
        for i in range(len(new_pairedbins)):
            if len(left) > 1:
                exons = left[0:len(left)-1] + new_pairedbins[i] + right
            else:
                exons = new_pairedbins[i] + right
            # Bin is added to the bin list
            new_bin = BinT(exons=exons, count=count)
            if (len(new_bin.exons) > 2):
                all_pairedbins.append(new_bin)
    return all_pairedbins



# GROUP PAIRS FUNCTION
"""
The group pairs function groups pairs with the identical first bin
"""
def group_pairs(pairedbins):
    group_dict = {} # The bins are stored in a dictionary
    for pairedbin in pairedbins:
        left = pairedbin.leftExons
        right = pairedbin.rightExons
        invalid = False
        
        if left[0] > right[0]:
            x = left
            left = right
            right = x
        j = 0
        for i in range(0,len(left)):
            if left[i] == right[j]:
                j += 1
                if len(right) > j+1:
                    if i!=(len(left)-1) and left[i+1] != right[j]:
                        invalid = True
                        break  
                else:
                    break
        if invalid == True:
            continue
        if len(right) == 0 or len(left) == 0:
            continue
        if right[0] < left[len(left)-1] and right[0] not in left:
            continue
        if len(right) == 1 and len(left) == 1 and left[0] + 1 == right[0]:
            continue
        if len(right) == 1 and right[0] in left:
            continue
        
        left = tuple(left) # Tuples as keys
        if left in group_dict.keys():
            right_exons = group_dict[left]
            if right not in right_exons:
                right_exons.append(right)
            group_dict[left] = right_exons
        else:
            right_exons = []
            right_exons.append(right)
            group_dict[left] = right_exons
    
    return group_dict

# FILTER_TRANSCRIPTS FUNCTION
"""
The filter_transcripts function removes all transcripts which match the common bin but none of the second partners (for each group)
"""
def filter_transcripts(transcripts,group_dict):
    for left_exons in group_dict.keys(): # Take each group
        i = 0
        while i < len(transcripts): # Iterate through the list of transcripts
            transcript = transcripts[i]
            validpath = True # Boolean checking if transcript is valid
            if left_exons[0] in transcript:
                j = transcript.index(left_exons[0]) # Get index where common bin starts
                if j+len(left_exons) <= len(transcript): # If this statement is false, the common bin is not part of the trancript
                    k = 1
                    valid = True # Boolean checking if the common bin is part of the transcript
                    while k < len(left_exons):
                        if transcript[j+k] != left_exons[k]:
                            valid = False 
                            break
                        k += 1
                    
                    if valid == True: # If the common bin is part of the transcript, let's look for the second partners
                        validpath = False
                        for right_exons in group_dict[left_exons]: # Iterate through second partners
                            if right_exons[0] in transcript:
                                start_index = transcript.index(right_exons[0])
                                if start_index+len(right_exons) <= len(transcript):
                                    m = 1
                                    valid2 = True
                                    while m < len(right_exons):
                                        if transcript[start_index+m] != right_exons[m]:
                                            valid2 = False
                                            break
                                        m += 1
                                    if valid2 == True:
                                        validpath = True
                                        break
                        if validpath == False: # If common bin was found in the transcript but none of the second partners, the transcript is removed
                            transcripts.pop(i)
                            i = i-1
            i+=1                  
    return transcripts
