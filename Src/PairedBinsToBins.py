#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple

from sqlalchemy import false
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript
from PathEnumeration import activeBinPathEnumeration3

def getBins(Bins:list):
    binList = []
    for bin in Bins:
        binList.append(bin[0])
    return binList

# Function to check, whether there are bins with incorrect order (e. g. left 1,3,2 or right 4,3,5)
def checkPairedBins(bin):
    for i in range (0, len(bin.leftExons)-2): 
        if bin.leftExons[i] > bin.leftExons[i+1]:
            return False 
    for i in range (0, len(bin.rightExons)-2): 
        if bin.rightExons[i] > bin.rightExons[i+1]: 
            return False            
    return True

def appendBinsLong(path_dict:dict, bin, fullLengthBoolean, Bins:list, allNewBins:list, BinT):
    for value in path_dict.values():
        if len(value)>=2:
            if fullLengthBoolean == True:
                newBin = bin.leftExons + value[1:] + bin.rightExons
            else:
                newBin = bin.leftExons + value[1:len(value)-2] + bin.rightExons
        else:
            newBin = bin.leftExons + bin.rightExons
        if appendBinsShort(Bins, newBin) == True:
            Bins.append(BinT(exons=newBin, count={}))
        if newBin not in allNewBins:
            allNewBins.append(newBin)

def appendBinsShort(Bins:list, exonBin:list):
    for bin1 in Bins: 
        if exonBin == bin1.exons:
            return False
    return True
       

def fromPairedBinsToBins(pairedBins, Bins, graph, Exons:list):
    BinT = namedtuple('BinT', 'exons count')
    allNewBins = []
    for bin in pairedBins:
        if checkPairedBins (bin) == False:
            continue                                                                                                              # Read Bins from PairedBins
        last_left_exon = bin.leftExons[len(bin.leftExons)-1]                                                                            # Define LastExon of left Bins and 
        first_right_exon = bin.rightExons[0]                                                                                            # FirstExon of RightBins
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
#   V   # 13. left 1,2,3    right 1
        
        while (first_right_exon < last_left_exon) and len(bin.rightExons)>0:                                                            # As long as first_right_exon is smaller than last_left_exon                                                                                                                              
            # Trim cases 3, 5, 6, 7, 8, 9, 10, 11, 12:
            if first_right_exon < bin.leftExons[0]:                                                                                     # if first_right_exon is smaller than the first exon of the left partner
                exonBin.append(first_right_exon)                                                                                        # Add it the list exonBin, containing exons smaller than the first left exon
                bin.rightExons.remove(bin.rightExons[0])                                                                                # remove from rightExons
                if len(bin.rightExons)>0:                                                                                               # if the rightExons list is not empty now
                    first_right_exon = bin.rightExons[0]                                                                                # assign the new first right exon to its variable 
                else:                                                                                                                   # otherweise
                    first_right_exon = -1                                                                                               # assign -1 to first_right_exon      
            
            # Discard case 4, 10 and 11 here:
            else:
                incorrectPairedBinBoolean = True                                                                                        # define Boolean to report incorrect PairedBins
                whileLoopBreakBoolean = False                                                                                           # define Boolean to exit While-loop, if necessary
                
                # Check if leftExons contain number of first_right_exon is
                                                                                                             
                for i in range (0, len(bin.leftExons)):                                                                                 
                    if first_right_exon == bin.leftExons[i]:                                                                            # Match: (e. g. case 11)
                        # Filter now case 11: left 1,2,5    right 0,1,2,8 (the hardest)                                                 
                        if i<len(bin.leftExons)-1 and len(bin.rightExons)>1:                                                            # check, whether there are additional elements in both, leftExons and rightExons
                            if bin.leftExons[i+1] != bin.rightExons[1]:                                                                 # successors miss-match
                                whileLoopBreakBoolean = True                                                                            # report, to exit whileLoop
                                incorrectPairedBinBoolean = True                                                                        # set report it via incorrectpairedBinBoolean
                            else:                                                                                                       # Successors match
                                incorrectPairedBinBoolean = False                                                                       # Therefore no incorrect PairedBin
                                bin.rightExons.remove(bin.rightExons[0])                                                                # Currently remove currently checked Exon
                                first_right_exon = bin.rightExons[0]                                                                    # Assign new first right exon to its variable
                            break                                                                                                       # Exit for-loop (since whileLoopBreakBoolean was initialized with False -> while-loop will not be exited)
                        else:                                                                                                           # no further list elements on at least one site 
                            incorrectPairedBinBoolean = False                                                                           # Therefore no incorrect PairedBin
                            whileLoopBreakBoolean = True                                                                                # Break whileLoop, because either: 
                                                                                                                                            # a. no further element in leftExons -> next element in right_exons > last_left_exon
                                                                                                                                            # b. no further element in rightExons -> no rightExon to check anymore  
                            bin.rightExons.remove(bin.rightExons[0])                                                                    # remove currently checked exon
                            if len(bin.rightExons)>0:                                                                                   # If rightExon list is not empty yet, 
                                first_right_exon = bin.rightExons[0]                                                                         # assign new first right exon to its variable
                            else:                                                                                                       # otherwise set first_right_exon to -1
                                first_right_exon = -1 
                            break                                                                                                       # Exit For-Loop
                if incorrectPairedBinBoolean == True or whileLoopBreakBoolean == True:                                                  # Exit while-loop if it's an incorrect PairedBin or While-LoopBreakBoolean has been set to True
                    break                                                                                                               
        if incorrectPairedBinBoolean == True:                                                                                           # If an incorrect PairedBin has been found, discard it 
            continue                                                                                                                    # Check next bin
                                                                                                                            
        # Cases that need to be further addressed: 0,1,2,3,5

        if last_left_exon == first_right_exon:                                                                                          # If the lastLeftExon is equal to the last rightExon we have to options
            if len(bin.leftExons)>1:                                                                                                        # a. if there's more than one element in leftExons: 
                bin.leftExons.remove(bin.leftExons[len(bin.leftExons)-1])                                                                        # remove last exon of leftExons and 
                last_left_exon = bin.leftExons[len(bin.leftExons)-1]                                                                             # assign new last left exon to its corresponding variable 
            elif len(bin.rightExons)>1:                                                                                                     # b. is there otherwise more than one element in rightExons and only one in leftExons 
                bin.rightExons.remove(bin.rightExons[0])                                                                                        # remove first right exon  
                first_right_exon = bin.rightExons[0]                                                                                            # and assign the new rightExon to its corresponding variable
            exonBin = exonBin + bin.leftExons + bin.rightExons                                                                              # exonBin might be empty or not (maybe we've already added exons from rightExons,
                                                                                                                                            # which were smaller than the first left exon), anyway, add remaining leftExons
                                                                                                                                            # and remaining rightExons   
            # Discard case 0 and 1 here:                                                                                                    
            if len(exonBin)>2:                                                                                                          # If exonBin contains more than two elements
                newBinBoolean = appendBinsShort (Bins, exonBin)                                                                                   # Check if Bins contain ExonBin already           
                if newBinBoolean == True:                                                                                                   # if this is a real Multi-Bin and it's not already in the Bins 
                    Bins.append(BinT(exons=exonBin, count={}))                                                                                  # add it to Bins
                if exonBin not in allNewBins:                                                                                               # if it´s not in the allNewBins list -> add it
                    allNewBins.append(exonBin)                                                                                               
            continue                                                                                                                    # check next Bin
        
        # Cases that have not or only partially been addressed yet: -1, 6, 7, 8, 12, 13
            # Case -1
        if len(bin.leftExons)==1 and len(bin.rightExons) == 1 and bin.leftExons[0] + 1 == bin.rightExons[0]:                            # If we have case -1 (e. g. left 0, right 1)-> nothing to do
            continue                                                                                                                    # Check next bin
        # Case 6, 7 and 8:
        if len(bin.rightExons) == 0:                                                                                                    # If rightExons were emptied previously
            #  6 and 7
            if len(exonBin)==0:                                                                                                         # and no Exons were added to exonBin
                if len(bin.leftExons)>2:                                                                                                # bin.leftExons is the remaining part of PairedBins -> if this is a real MultiBin
                    if appendBinsShort(Bins, bin.leftExons) == True:                                                                        # check whether Bins already contain it's alrdy contained in Bins
                        Bins.append(BinT(exons=bin.leftExons, count={}))                                                                    # add it, if necessary
                    if bin.leftExons not in allNewBins:                                                                                 # add it to allNewBins, if it's not in there alrdy 
                        allNewBins.append(bin.leftExons)
                continue                                                                                                                # check next bin
            else:                                                                                                                       # Otherwise exonBin is not empty
                if exonBin[len(exonBin)-1] == bin.leftExons[0]:                                                                             # if last exon of exonBin and first exon of leftExons are the same
                    exonBin = exonBin + bin.leftExons[1:]                                                                                       # add all exons from leftExons, but not the first 
                    if appendBinsShort(Bins, exonBin) == True:                                                                                  # Check, whether Bins contains this pairedBin already
                        Bins.append(BinT(exons=exonBin.leftExons, count={}))                                                                        # add it, if its new
                    if exonBin not in allNewBins:                                                                                               # add it also to allNewBins, if it's not already in there
                        allNewBins.append(exonBin)
            # Case 8
                else:                                                                                                                       # Otherwise last exon of exonBin and first exon of leftExons are not the same 
                    startNodeLastLeftExon = exonBin[len(exonBin)-1] + 2                                                                         # define startNode of last_left_exon
                    startNodeFirstRightExon = bin.leftExons[0] + 2                                                                              # define rightnode of last_right_exon
                    full_path_dict = activeBinPathEnumeration3(str(startNodeFirstRightExon), str(startNodeLastLeftExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins) # employ path enumeartion using the multiBinconstraints
                    if full_path_dict != None:                                                                                                  # if there was at least one path found
                        appendBinsLong(full_path_dict, bin, False, Bins, allNewBins, BinT)                                                          # append it to the dictionary
        
        if last_left_exon < first_right_exon:                                                                                               # last case -> exons left in rightExons and last_left_exon<first_right_exon  
            startNodeLastLeftExon = last_left_exon + 2                                                                                          # define startnode of last_left_exon
            startNodeFirstRightExon = first_right_exon +2                                                                                       # Get startnode of firstRightExon                                                                                
            full_path_dict = activeBinPathEnumeration3(str(startNodeFirstRightExon), str(startNodeLastLeftExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins) # employ pathEnumeration with MultiBinConstraint
            if full_path_dict != None:                                                                                                      # If at least one path between real_last_left_exon and real_first_right_exon has been found
                appendBinsLong(full_path_dict, bin, True, Bins, allNewBins, BinT)                                                           # append it
    return Bins, allNewBins