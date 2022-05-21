#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript

def getBins(Bins:list):
    binList = []
    for bin in Bins:
        binList.append(bin[0])
    return binList

def checkPairedBins(PairedBins):
    binList = []
    for bin in PairedBins:
        incorrectPairedBinBoolean = False
        for i in range (0, len(bin.leftExons)-1): 
            if incorrectPairedBinBoolean == True:
                break
            for j in range (i+1, len(bin.leftExons)):
                if bin.leftExons[i] > bin.leftExons[j]: 
                    incorrectPairedBinBoolean = True
                    break
                
        for i in range (0, len(bin.rightExons)-1): 
            if incorrectPairedBinBoolean == True:
                break
            for j in range (i+1, len(bin.rightExons)):
                if bin.rightExons[i] > bin.rightExons[j]: 
                    incorrectPairedBinBoolean = True
        if incorrectPairedBinBoolean == False:
            binList.append(bin)
    return binList
        

def fromPairedBinsToBins(pairedBins, Bins, graph, Exons:list):
    BinT = namedtuple('BinT', 'exons count')
    allNewBins = []
    # Read Bins from PairedBins
    for bin in pairedBins:
        #If this the first ListElement of PairedBins, create a List, to save all Bins, that have already been added
        # Define LastExon of left Bins and FirstExon of RightBins
        real_last_left_exon = bin.leftExons[len(bin.leftExons)-1]
        real_first_right_exon = bin.rightExons[0]
        # Make a working-copy of these Exons 
        last_left_exon = real_last_left_exon
        first_right_exon = real_first_right_exon
        incorrectPairedBinBoolean = False                                                                                               # Define the incorrectPairedBinBoolean
        exonBin = []                                                                                                                    # Define new exonBin List for the Bins, which is later returned
        
        # Cases we need to adress:
        # -1 left 0         right 1
#   I   # 0. left 1         right 1
#   I   # 1. left 1,2       right 2   -> can be easily handled within the if last_left_exon == first_right_exon block and doesn't need path enumeration and can be discarded because it's not a real MultiBin 
#   V   # 2. left 1,2       right 2,5 -> can be easily handled within the if last_left_exon == first_right_exon block and there's no need for path-enumeration, but it can be added to Bins, because its a real MultiBin
#   V   # 3. left 1,2,5     right 2,5,8 -> 2 & 5 from rightExons need to be eliminated, but no path-enumeration is needed                            
#   I   # 4. left 1,2,5     right 3,5,8 -> needs to be discarded (because 3-5 is not compatible with 2-5)
#   V   # 5. left 1,2,5     right 1,2   -> needs to be trimmed down to 1,2,5 and added as MultiBin, if its not already in the Bins
#   V   # 6. left 1,2,5     right 0   -> needs to be transformed to 0,1,2,5
#   V   # 7. left 2,5       right 0,1 -> needs to be transformed to 0,1,2,5 -> no path enumeration needed
#   V   # 8. left 4,5       right 0,1 -> pathEnumeration needed from 1 to 4
#   V   # 9. left 2,5       right 0,1,2,5,8 -> combined case from 3, 5 and 6 -> no path enumeration needed
#   I   # 10. left 2,5      right 0,1,3,5,8 -> combined case from 4 and 6 -> needs to be discarded
#   I   # 11. left 1,2,5    right 0,1,2,8 -> the hardest -> needs to be discarded
#   V   # 12. left 1,2,3    right 5,7,8 -> the easiest -> actually path enumeration
        # First Step: Rearrange cases with mixed-up paired-bins
        # Case 5 to: left 1,2   right 1,2,5 
        # Case 6 to: left 0     right 1,2,5
        # Case 7 to: left 0,1   right 2,5
        # Case 8 to: left 0,1   right 4,5
        # Next Trimm Exons not needed
        while (first_right_exon < last_left_exon) and len(bin.rightExons)>0:                                                           # As long as first_right_exon is smaller than last_left_exon                                                                                                                              
            # Trim down the following cases:
            # 3.  left 1,2,5    right 2,5,8
            # 5.  left 1,2,5    right 1,2
            # 6.  left 1,2,5    right 0
            # 7.  left 2,5      right 0,1
            # 8.  left 4,5      right 0,1
            # 9.  left 2,5      right 0,1,2,5,8
            # 10. left 2,5      right 0,1,3,5,8
            # 11. left 1,2,5    right 0,1,2,8
            # 12. left 1,2,3    right 5,7,8 
            if first_right_exon < bin.leftExons[0]:                                                                                     # if first_right_exon is smaller than the first exon of the left partner
                print(bin)
                exonBin.append(first_right_exon)
                bin.rightExons.remove(bin.rightExons[0])                                                                                # remove this exon from rightExons
                if len(bin.rightExons)>0:                                                                                               # if the rightExons list is not empty
                    first_right_exon = bin.rightExons[0]                                                                                # assign the new first right exon to its variable 
                else:                                                                                                                   # otherweise
                    first_right_exon = -1                                                                                               # assign -1 to first_right_exon      
            # 4. left 1,2,5     right 3,5,8 -> discard it here
            # 10. left 2,5      right 0,1,3,5,8 -> discard it here
            else:                                                                                                                       # if first_right_exon is not smaller than the first exon of the left partner
                if first_right_exon not in bin.leftExons:                                                                               # check whether this exon is in bin.leftExons
                    incorrectPairedBinBoolean = True                                                                                    # if not, report by setting the incorrectPairedBinBoolean to true, 
                                                                                                                                        # which will be used later to discard this bin
                    break                                                                                                               # escape from while-loop
                else:                                                                                                                   # otherwise you found an exon in the left partner with the same value as the 
                    bin.rightExons.remove(bin.rightExons[0])                                                                            # currently check exon -> which remove this Exon and
                                                                                                                                                                                                   
                    if len(bin.rightExons)>0:                                                                                           # if rightExon list is not empty yet, 
                        first_right_exon = bin.rightExons[0]                                                                            # assign new first right exon to its variable
                    else:                                                                                                               # otherwise set first_right_exon to -1
                        first_right_exon = -1    

        if incorrectPairedBinBoolean == True:                                                                                           # If an incorrect PairedBin has been found, discard it 
            continue                                                                                                                    # Check next bin
        # Cases that need further to be addressed in the upcoming part:
            # 0. left 1         right 1
            # 1. left 1,2       right 2
            # 2. left 1,2       right 2,5
            # 3. left 1,2,5     right 2,5,8
            # 5. left 1,2,5     right 1,2

        if last_left_exon == first_right_exon:                                                                                          # If the lastLeftExon is equal to the last rightExon we have to options
            if len(bin.leftExons)>1:                                                                                                    # if there's more than one element in leftExons: 
                bin.leftExons.remove(bin.leftExons[len(bin.leftExons)-1])                                                               # remove last exon of leftExons and 
                last_left_exon = bin.leftExons[len(bin.leftExons)-1]                                                                    # assign new last left exon to its corresponding variable 
            elif len(bin.rightExons)>1:                                                                                                 # is there otherwise more than one element in rightExons and only one in leftExons 
                bin.rightExons.remove(bin.rightExons[0])                                                                                # remove first right exon  
                first_right_exon = bin.rightExons[0]                                                                                    # and assign the new rightExon to its corresponding variable
            
            exonBin = exonBin + bin.leftExons + bin.rightExons                                                                           # exonBin might be empty or not (maybe we've already added exons from rightExons,
                                                                                                                                        # which were smaller than the first left exon), anyway, add remaining leftExons
            # Discard 0. and 1. here                                                                                                    # and remaining rightExons  
            if len(exonBin)>2:
                newBinBoolean = False
                for bin1 in Bins: 
                    if exonBin == bin1.exons:
                        newBinBoolean=False
                        break
                if newBinBoolean == True:                                                                                              # if this is a real Multi-Bin and it's not already in the Bins 
                    Bins.append(BinT(exons=exonBin, count={}))
                if exonBin not in allNewBins:
                    allNewBins.append(exonBin)                                                                                              # add it to Bins 
            continue                                                                                                                    # check next Bin
        
        # Cases that have not or only partially been addressed yet:
            # -1. left 0        right 1
            # 6.  left 1,2,5    right 0   -> needs to be transformed to 0,1,2,5
            # 7.  left 2,5      right 0,1 -> needs to be transformed to 0,1,2,5 -> no path enumeration needed
            # 8.  left 4,5      right 0,1 -> pathEnumeration needed from 1 to 4
            # 11. left 1,2,5    right 0,1,2,8 -> the hardest
            # 12. left 1,2,3    right 5,7,8 -> the easiest 
        # 13. left 1,2,3    right 1     
        # Case -1
        if len(bin.leftExons)==1 and len(bin.rightExons) == 1 and bin.leftExons[0] + 1 == bin.rightExons[0]:                               # If we have case -1 (e. g. left 0, right 1)-> nothing to do
            continue
        # Case 6, 7 and 8:
        if len(bin.rightExons) == 0:
            #  6 and 7
            if len(exonBin)==0:
                if len(bin.leftExons)>2:
                    newBinBoolean = True
                    for bin1 in Bins: 
                        if bin.leftExons==bin1.exons:
                            newBinBoolean = False
                            break
                    if newBinBoolean == True:
                        Bins.append(BinT(exons=bin.leftExons, count={}))
                    if bin.leftExons not in allNewBins:
                        allNewBins.append(bin.leftExons)
                continue
            else: 
                newBinBoolean = True
                if exonBin[len(exonBin)-1] == bin.leftExons[0] + 1:
                    startNodeLastLeftExon = exonBin[len(exonBin)-1] + 2  
                    startNodeFirstRightExon = bin.leftExons[0] + 2
                    full_path_dict = pairedBinToBins2(str(startNodeLastLeftExon), str(startNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins)
                    if full_path_dict != None:
                        for value in full_path_dict.values():
                            if len(value)>2:
                                newBin = bin.leftExons + value[1:len(value)-2] + bin.rightExons
                            else:
                                newBin = bin.leftExons + bin.rightExons
                        # If this Bin is not already in Bins
                            for bin1 in Bins:
                                if newBin==bin1.exons:
                                    newBinBoolean = False
                                    break
                            if newBinBoolean == True:
                                Bins.append(BinT(exons=newBin, count={}))
                            if newBin not in allNewBins:
                                allNewBins.append(newBin)
                    #exonBin = exonBin + bin.leftExons
                    # PathEnumeration?
                    # for bin1 in Bins:
                    #     if exonBin == bin1.exons:
                    #         newBinBoolean = False
                    #         break
                    # if newBinBoolean == True:
                    #     Bins.append(BinT(exons=exonBin.leftExons, count={}))
                    # if exonBin not in allNewBins:
                    #     allNewBins.append(exonBin)
                elif exonBin[len(exonBin)-1] == bin.leftExons[0]:
                    exonBin = exonBin + bin.leftExons[1:]
                    for bin1 in Bins:
                        if exonBin == bin1.exons:
                            newBinBoolean = False
                            break
                    if newBinBoolean == True:
                        Bins.append(BinT(exons=exonBin.leftExons, count={}))
                    if exonBin not in allNewBins:
                        allNewBins.append(exonBin)
            # Case 8
                else:
                    startNodeLastLeftExon = exonBin[len(exonBin)-1] + 2  
                    startNodeFirstRightExon = bin.leftExons[0] + 2
                    full_path_dict = pairedBinToBins2(str(startNodeLastLeftExon), str(startNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins)
                    if full_path_dict != None:
                        for value in full_path_dict.values():
                            if len(value)>2:
                                newBin = bin.leftExons + value[1:len(value)-2] + bin.rightExons
                            else:
                                newBin = bin.leftExons + bin.rightExons
                        # If this Bin is not already in Bins
                            for bin1 in Bins:
                                if newBin==bin1.exons:
                                    newBinBoolean = False
                                    break
                            if newBinBoolean == True:
                                Bins.append(BinT(exons=newBin, count={}))
                            if newBin not in allNewBins:
                                allNewBins.append(newBin)
        
        if last_left_exon<first_right_exon:                                                                                             # Check if this particular exonBin has allred been added to the list 
            startNodeLastLeftExon = last_left_exon + 2
            startNodeFirstRightExon = first_right_exon +2                                                                               # Get endNode of current firstRightExon                                                                                
            full_path_dict = pairedBinToBins2(str(startNodeLastLeftExon), str(startNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins)
            #If at least one path between real_last_left_exon and real_first_right_exon has been found
            if full_path_dict != None:
                for value in full_path_dict.values():
                    newBinBoolean = True
                    if len(value)>=2:
                        newBin = bin.leftExons + value[1:] + bin.rightExons
                    else:
                        newBin = bin.leftExons + bin.rightExons
                    # If this Bin is not already in Bins
                    for bin1 in Bins:
                        if newBin==bin1.exons:
                            newBinBoolean = False
                            break
                    if newBinBoolean == True:
                        Bins.append(BinT(exons=newBin, count={}))
                    if newBin not in allNewBins:
                        allNewBins.append(newBin)
    return Bins, allNewBins

def pairedBinToBins2(v:str, endNode:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v not in list(graph.nodes()):
        return
    if v == endNode:                                                                                                        # If current node is drain
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
                # """ else:
                #     if bin1.exons[0] == activeExonNumber:                                                           
                #         newActiveBins.append(bin1)                                                                                            
                #  """                                                                                                        
        else:                                                                                                           # If EdgeType is no SpliceJunctions
            newActiveBins=activeBins                                                                                    # all formerly activeBins are newActiveBins
        pairedBinToBins2(u, endNode, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)                       # call ActiveBinPathEnumeration for the current node
    return allpaths

def PairedBinToBinEnumeration (v:str, endNode:str, pfad:list, allpaths:dict, path_number:list, graph):  
    if v not in list(graph.nodes()):
        return
    if v == endNode:
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
        path_number[0] = path_number[0] + 1
        return
    else:
        for u in graph.adj[str(v)].keys():
            PairedBinToBinEnumeration(u, str(endNode), pfad+[u], allpaths, path_number, graph)
    return allpaths
