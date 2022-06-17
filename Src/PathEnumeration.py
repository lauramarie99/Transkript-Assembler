#!/usr/bin/env python

import sys, ast, os
from collections import namedtuple

from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript
   
# Get new Bins without MultiBinRestriction
def getNewActiveBinsEasily (newActiveBins:list, Bins:list, current_Exon:str):
    for bin in Bins:    
        if (bin.exons[0])==current_Exon:
            newActiveBins.append(bin)

# Get new Bins with MultiBinRestriction
def getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, act_Exon, SuccInMultiBinsList):
    for bin in Bins:
        if bin.exons[0]==act_Exon:
            # If it's a Bin containing only 1 Exon add it to newActiveBins
            if len(bin.exons) == 1: 
                newActiveBins.append(bin)
            else: 
                # If it's a bin containing at least two exons 
                if bin.exons[1] in SuccInMultiBinsList:       
                    # If the exon following the current exon in this Bin is represented by a currently active MultiBin add the Bin to the newActive Bin List       
                    newActiveBins.append(bin)

#Get only Bins with length>2
def getMultiBins(Bins:list):
    multiBins = []
    for bin in Bins:
        if len(bin.exons)>2:
            multiBins.append(bin)
    return multiBins

# Full Path Enumeration
def fullPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, graph):  
    if v == '1':                                                                                                        # If current is drain, return current path
        allpaths[path_number[0]] = pfad                                                                                 # Write current path transcribed into a transcript in to the 
                                                                                                                        # allPath-Dictionary
        path_number[0] = path_number[0] + 1                                                                             # Increase PathNumber by one
        return
    else:
        for u in graph.adj[str(v)].keys():                                                                              #For all successors of the current node repeat this procedure
            fullPathEnumeration(u, pfad+[u], allpaths, path_number, graph)
    return allpaths

#Enumerate Paths with for-Loop later-on -> Problem: more if-Statements required to control all the exceptions
# All Bins are written into ActiveBins.
def activeBinPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    # If current node ist drain, return current path
    if v == '1':                                                                                                        # If current node is drain
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)                                                  # Add current path transcribed into a transcript to the allpath dictionary
        path_number[0] = path_number[0] + 1                                                                             # Increase Pathcounter
        return
    kantenTyp = None                                                                                                    # Define Variable for Kantentyp
    if len(pfad)>1:                                                                                                     # If pathlength is great than 1 (node is not source), retrieve edgetype
        kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']                                                      
    if kantenTyp == 'SpliceJunction':                                                                                   # Only if edgeType is a SpliceJunction
        activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']                                              # Define beginning of the edge as previousExonNumber
        previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']                                          # Define end of the edge as ActiveExonNumbger
        newActiveBins = []                                                                                              # Make a new empty List for the newActiveBins
        SuccInMultiBinsList = []                                                                                        # Make a new List for Successors in MultiBins
        compatibleBinBoolean = False                                                                                    # Define a new Boolean to check for compatible Bins in all ActiveBins
        potentiallyConflictingBinBoolean = False                                                                        # Define a secondBoolean, to check whether there are conflicting Booleans with 
        for bin in activeBins:
            for j in range(1, len(bin.exons)):                                                                          # Check for compatible Bins                                                     
                if bin.exons[j] == activeExonNumber and bin.exons[j-1] == previousExonNumber:                           # If you found at least one compatible Bin 
                    if j-1>0:
                        compatibleBinBoolean = True                                                                     # Report, that there are compatibleBins
                    if j<len(bin.exons)-1:                                                                              # If the end of the current edge is not the last Exon of the current Bin
                        newActiveBins.append(bin)                                                                       # Add this Bin to the newActiveBins                                                                                                                                
                        SuccInMultiBinsList.append(bin.exons[j+1])                                                      # Add Successors of the currently active Exon (end ofo the current edge) to a succesorInMultBinList 
                elif bin.exons[j-1] == previousExonNumber and j-1>0:                                                    # If you find potentially conflicting bins, (e. g. current edge 2-4 and the currently 
                                                                                                                        # checked bin is 1-2-3)
                    potentiallyConflictingBinBoolean = True                                                             # Report that there are potentially conflicting Exons
                                                                                                                                            
        if compatibleBinBoolean == False and potentiallyConflictingBinBoolean == True:                                  # If there are no compatibleExons, whose firstExon is not the beginning of the current edge, but there 
            return                                                                                                      # are conflicting Exons, stop exeecuting this path.
        if len(SuccInMultiBinsList) == 0:                                                                               # If you SuccessorInMultiBinlist is empty, GetNewActiveBins without restrictions
            getNewActiveBinsEasily(newActiveBins, Bins, activeExonNumber)   
        else:                                                                                                           # Else: use a more complex way to retrieve the newActiveBins
            getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, activeExonNumber, SuccInMultiBinsList)        
    # If edgeType is not a SpliceJunction
    else:
        newActiveBins=activeBins                                                                                        # If its not a SpliceJunction, all activeBins are maintained
    for u in graph.adj[str(v)].keys():                                                                                  # Iterate through all successors of the current Node and 
        activeBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)                        # Call ActivePathEnumeration
    return allpaths

#Enumerate Paths with ForLoop initially
#All Bins/MultiBins are transferred to the Variable activeBins
def activeBinPathEnumeration2(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v == '1':                                                                                                            # If current node ist drain
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)                                                      # Add current path transcribed into a transcript to the allpath dictionary
        path_number[0] = path_number[0] + 1                                                                                 # Increase Pathcounter
        return
    for u in graph.adj[str(v)].keys():                                                                                                 
        if graph.edges[v, u]['type'] == 'SpliceJunction':                                                                   # If current Edge was a SpliceJunction
            activeExonNumber=graph.edges[v, u]['endExon']                                                                   # Define beginning of the edge as previousExonNumber 
            previousExonNumber=graph.edges[v,u]['startExon']                                                                # Define end of the edge as ActiveExonNumbger
            newActiveBins = []                                                                                              # Make a new empty List for the newActiveBins
            compatibleBinBoolean = False                                                                                    # Define a new Boolean to check for compatible Bins in all ActiveBins
            potentiallyConflictingBinBoolean = False                                                                        # Define a secondBoolean, to check whether there are conflicting Booleans with 
                                                                                                                            # the ones compatible to the current path (e. g. CurrentPath: 1->2->4, currentActiveBoolean 2,4,7;
                                                                                                                            # furtherActiveBoolean is 1,2,3, but there#s no 1-2-4) 
            for bin in activeBins:                                                                                          # Check active Bins now for compatible Bins
                for j in range(0, len(bin.exons)-1):                                                                        # Check all positions in the ActiveBin
                    if bin.exons[j+1] == activeExonNumber and bin.exons[j] == previousExonNumber:                           # If the current Edge is represented in this Boolean:
                                                                                                                            
                        if j!=0:                                                                                            # Make sure, that the beginning of this edge is not the beginning of this Bin (e.g. 2-4-7). 
                                                                                                                            # Otherwise there might be conflicting bins (such as 1,2,3 but without 1,2,4 available).
                                                                                                                            # This could result in abrogation of following this path. If there are no conflicting booleans also OK.                                                                                                     
                            compatibleBinBoolean = True                                                                     # If it's not (e. g. current edge is 2-4 and currently checked active Bin is 1-2-4, this path is valid!!!), 
                                                                                                                            # because other Bins, such as 1-2-3 are not conflicting in this case. 
                        if j!=len(bin.exons)-2:                                                                             # Make sure current ActiveExon (end of the edge) is not the Endexon in the currentlyChecked bin 
                                                                                                                            # (otherwise it would no be required anymore).
                            newActiveBins.append(bin)                                                                       # Anyway, write this Bin into newActiveBins, in case this path is followed
                    elif (bin.exons[j] == previousExonNumber) and j!=0:                                                     # Check for conflicting Bins (e. g. current edge 2-4 with Bin 2-4-7), but 1-2-3 is in ActiveBins, then 1-2-4 
                                                                                                                            # would be required as well, otherwise following the current path should be stopped.   
                        potentiallyConflictingBinBoolean = True                                                             # Anyway, if there is another Bin containing the StartExon of the currentEdge (previousExonNumber), 
                                                                                                                            # report, that there are potentially conflicting Bins 

            if compatibleBinBoolean == False and potentiallyConflictingBinBoolean==True:                                    # If there are Bins that are in conflict with the currently followed path, cancel following this path
                continue                                                                                                    # Check next Edge
            for bin1 in Bins:                                                                                               # Check for new active Bins
                if bin1.exons[0] == activeExonNumber:                                                                       # If there's a Bin in all Multibins, whose starting Exon ist the end of the current Edge, 
                    newActiveBins.append(bin1)                                                                              # Add this bin to the MultiBinList
        else:                                                                                                               # If the current Edge is not a SpliceJunction                            
            newActiveBins=activeBins                                                                                        # All ActiveBins are newActiveBins.
        activeBinPathEnumeration2(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)                           # Call activeBinPathEnumeration2 for the current Node 
    return allpaths                                                                                                         # In the end, return the dictionary of all GeneX pathdictionaries with enumerated paths 0,1,2,....

# Third possibility to Enumerate, actually the easiest and smoothest one
# Main important difference to ActiveBinPathEnumeration2:
    # ActiveBins are only considered as such if the end of the current edge (ActiveExonNumber) is not the first Exon in the currently cheked Bin
    # This ensures the following:
        # 1. Only Exons are added, newly, to the newActiveBins (already active and valid exons are kept anyway), whose second exon is 
        #    the end of the current edge. Therefore, no check is necessary anymore, whether the beginning of the upcoming edges 
        #    (previousExonNumber), are the firstExon in the ActiveBins checked for all following nodes. 
        # 2. The check for conflicting bins is then obsolete, because, if there were (e. g. 2->4 (for Exon 2 -> 1,2,3 was added, 2-4-7 
        #    was not added yet, because 2 was not a middle-exon) and no other activeBins were available (e.g. 1-2-4 is missing) no 
        #    compatible would be found. Otherwise, if there were other compatible Bins (such as 1-2-4) those would have been added in 
        #    the previous iteration
    # Making sure, that not all paths are discared, an empty ActiveBinList is written into ActiveBins during fist call of 
    # activeBinpathEnumeration3     
    #  

def activeBinPathEnumeration3(endnode:str, v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v not in list(graph.nodes()):
        return
    if v == endnode:                                                                                                    # If current node is drain
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)                                                  # Add current path transcribed into a transcript to the allpath dictionary
        path_number[0] = path_number[0] + 1                                                                             # Increase Pathcounter
        return
    for u in graph.adj[str(v)]:                                                                                         # For all successors of the current node
        if graph.edges[v, u]['type'] == 'SpliceJunction':                                                               # If edgeType is a SpliceJunction
            newActiveBins = []                                                                                          # Create an empty List for new Active Bins
            activeExonNumber=graph.edges[v, u]['endExon']                                                               # Define beginning of the edge as previousExonNumber
            previousExonNumber=graph.edges[v,u]['startExon']                                                            # Define end of the edge as ActiveExonNumbger
            compatibleBinListBoolean = False                                                                            # Define a new Boolean to check for compatible Bins in all ActiveBins
            for bin in activeBins:                                                                                      # For all ActiveBins
                for j in range(0, len(bin.exons)-1):                                                                    # For alle Positions in ActiveBins except the last
                    if bin.exons[j+1] == activeExonNumber and bin.exons[j] == previousExonNumber:                       # If you found at least one compatible Bin
                        compatibleBinListBoolean = True                                                                 # Report finding of an compatible bin (remember this Bin's first Exon is not the beginning of the current edge)
                                                                                                                        # e.g. curreent edge is 2-4 and current path is 1-2-4 (1,2,4 is in ActiveBinList, 1-2-3 as well, but not 2,4,7 or 0-2-4)
                        if j<len(bin.exons)-2:                                                                          # If currently ActiveExon is not the end of the Bin
                            newActiveBins.append(bin)                                                                   # add it to newActiveBins
                            
            if compatibleBinListBoolean == False and len(activeBins)>0:                                                 # If there are no compatible Bins, but the ActiveBinList was not empty
                continue                                                                                                # Skip this path
            for bin1 in Bins:                                                                                           # Check all Bins for newActive Bins.
                if len(bin1.exons)>2 and bin1.exons[1] == activeExonNumber and bin1.exons[0] == previousExonNumber:     # newly added bins to NewActiveBins are defined as such:
                    newActiveBins.append(bin1)                                                                          # 1. first Exon of the Bin is the beginning of the currently checked edge (PreviousExonNumber)
                                                                                                                        # 2. second Exon of the Bin is the end of the currently checked edge (ActiveExonNumber)
                                                                                                                        # If these criteria are met, add the bin to newActiveBins                                                                                                    
        else:                                                                                                           # If EdgeType is no SpliceJunctions
            newActiveBins=activeBins                                                                                    # all formerly activeBins are newActiveBins
        activeBinPathEnumeration3(endnode, u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)              # call ActiveBinPathEnumeration for the current node
    return allpaths