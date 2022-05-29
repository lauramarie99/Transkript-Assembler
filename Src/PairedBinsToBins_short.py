#!/usr/bin/env python

from itertools import count
from re import X
import sys, ast, os
from xmlrpc.client import Boolean, boolean
import networkx as nx
from collections import namedtuple

from sqlalchemy import false
from sympy import continued_fraction
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

def appendBinsLong(path_dict:dict, bin, fullLengthBoolean, Bins:list, allNewBins:list, BinT, count, leftExons, rightExons):
    for value in path_dict.values():
        if len(value)>=2:
            newBin = leftExons + value[1:len(value)-1] + rightExons
        else:
            newBin = leftExons + rightExons
        if len(newBin)<=2:
            continue
        if appendBinsShort(Bins, newBin, allNewBins) == True:
            allNewBins.append(BinT(exons=newBin, count=count))

def appendBinsShort(Bins:list, exons:list, allNewBins):
    for bin1 in Bins: 
        if exons == bin1.exons:
            return False
    if len(allNewBins)>0:
        for bin2 in allNewBins:
            if exons==bin2.exons:
                return False
    return True
            
def fromPairedBinsToBinsShort(pairedBins, Bins, graph, Exons:list, multiBinRestriction):
    BinT = namedtuple('BinT', 'exons count')
    allNewBins = []
    for bin in pairedBins:                                                                                                                # Read Bins from PairedBins
        
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
#   I   # 14. left 1,4,3    right 5,6,7 -> false order
#   i   # 15. left 4,6      right 2,5   -> needs to be discarded

    #1. Eliminate Cases with potential order mistakes
        # e. g. Case 14 
        if checkPairedBins (bin) == False:
            print(bin)                                                                                              # Check potential order mistakes
            continue
        rightExons = bin.rightExons                                                                                 # Assign bin.rightExons to rightExons
        leftExons = bin.leftExons                                                                                   # Assign bin.leftExons to leftExons
        count = bin.count                                                                                           # Assign bin.count to leftExons

    # 2. Switch leftExons and rightExons, if firstExon of rightExons is smaller than firstExon of leftExons
        # Case 6: left 1,2,5    right 0         ->  left 0          right 1,2,5
        # Case 7: left 2,5      right 0,1       ->  left 0,1        right 2,5
        # Case 8: left 4,5      right 0,1       ->  left 0,1        right 4,5
        # Case 9: left 2,5      right 0,1,2,5,8 ->  left 0,1,2,5,8  right 2,5
        # Case 10:left 2,5      right 0,1,3,5,8 ->  left 0,1,3,5,8  right 2,5
        # Case 11:left 1,2,5    right 0,1,2,8   ->  left 0,1,2,8    right 1,2,5
        # Case 15:left 4,6      right 2,5       ->  left 2,5        right 4,6
         
        if rightExons [0] < leftExons[0]:
            switchExons = rightExons
            rightExons = leftExons
            leftExons = switchExons
    
    # 3. Eliminate cases, in which the first rightExon is not in leftExons, but smaller than last leftExons
        #          case 4    left 1,2,5      right 3,5,8
        #          case 15   left 2,5        right 4,6
        # switched case 10:  left 0,1,3,5,8  right 2,5

        if rightExons[0]<leftExons[len(leftExons)-1] and rightExons[0] not in leftExons:                                          
            continue

    # 4. Trim rightExons:    

        # Case 0    left 1          right 1         -> left 1           right []
        # Case 1    left 1,2        right 2         -> left 1,2         right []
        # Case 2    left 1,2        right 2,5       -> left 1,2         right 5
        # Case 3    left 1,2,5      right 2,5,8     -> left 1,2,5       right 8
        # Case 4.1  left 1,2,5      right 2,3,5,8   -> left 1,2,5       right 3,5,8 -> discarded
        # Case 5    left 1,2,5      right 1,2       -> left 1,2,5       right []
# switched Case 9   left 0,1,2,5,8  right 2,5       -> left 0,1,2,5,8   right []
# switched Case 11  left 0,1,2,8    right 1,2,5     -> left 0,1,2,8     right 5 -> discarded
        # Case 13   left 1,2,3      right 1         -> left 1,2,3       right []

        incorrectPairedBinBoolean = False
        for i in range (0, len(leftExons)):                                                                                    
            if rightExons[0] == leftExons[i]:                                                                                   # Match: (e. g. case 11)
                # Filter now switched  case 11:     right 0,1,2,8   left 1,2,5 (the hardest)
                #        and switched case 4.1:     right 1,2,5     left 2,3,5,8                                                      
                rightExons.pop(0)                                                                                               # Remove firstRightExon
                if len(rightExons)>0:                                                                                           # If there's a remaining Exon in rightExons
                    if i<len(leftExons)-1 and leftExons[i+1] != rightExons[0]:                                                  # check, whether there are additional elements in both, leftExons and rightExons                                                                                             # successors miss-match
                        incorrectPairedBinBoolean = True                                                                        # report it via incorrectpairedBinBoolean
                        break                                               
                else:                                                                                                           # If rightExons is empty exit for-loop
                    break                                                                                                       # Exit for-loop (since whileLoopBreakBoolean was initialized with False -> while-loop will not be exited)                                                                                                                                                    

        if incorrectPairedBinBoolean == True:                                                                                   # If an incorrect PairedBin has been found, discard it 
            continue
        
    # 5. Empty rightExons don't need pathEnumration
        #              Trimmed Case 0   left 1           right []
        #              Trimmed Case 1   left 1,2         right []
        #              Trimmed Case 5   left 1,2,5       right []
        # Switched and Trimmed Case 9   left 0,1,2,5,8   right []
        #              Trimmed Case 13  left 1,2,3       right []
        
        if len(rightExons) == 0:
            # Only if it's a real MultiBin -> Discard Case 0,1
            if len(leftExons)>2:
                if appendBinsShort(Bins, leftExons, allNewBins) == True:                                                        # check whether Bins already contain it's alrdy contained in Bins
                    allNewBins.append(BinT(exons=leftExons, count=count))                                                       # add it, if necessary
            continue                                                                                                            # Check next bin
    
    # 6. Deal with all cases with no empty rightExons
        #           Case -1     left 0          right 1
        # Trimmed   Case 2      left 1,2        right 5
        # Trimmed   Case 3      left 1,2,5      right 8                                                                                                            
        # Switched  Case 6      left 0          right 1,2,5
        # Switched  Case 7      left 0,1        right 2,5
        # Switched  Case 8      left 0,1        right 4,5
        #           Case 12     left 1,2,3      right 5,7,8
        else:                                                                                                                   # Check Bins where len(rightExons)>1
            if multiBinRestriction == True:                                                                                     # If MultiBinRestriction is on -> check that resulting length is >2                               
            # 6.1 Get rid of case -1
                if len(leftExons)==1 and len(rightExons) == 1 and leftExons[0] + 1 == rightExons[0]:                            # If we have case -1 (e. g. left 0, right 1)-> nothing to do
                    continue

            # 6.2 Now enumerate all residual cases   
                # Trimmed Case 2:       left 1,2         right 5
                # Trimmed Case 3:       left 1,2,5       right 8
                # Switched Case 6:      left 0           right 1,2,5
                # Switched  Case 7      left 0,1        right 2,5
                # Switched  Case 8      left 0,1        right 4,5
                #           Case 12     left 1,2,3      right 5,7,8
                
            startNodeLastLeftExon = -1
            endNodeFirstRightExon = -1
            
            # 6.2.1 extract beginning end ending nodes of left and rightExon safely  
            if leftExons[len(leftExons)-1] < rightExons[0]:                                                                     
                for edgeKey, edgeValue in graph.edges.items():                                                                  # Iterate over all edges in the Graph
                    if edgeValue['type'] == 'Exon':                                                                             # If you found an Exon 
                        if edgeValue['exon'] == leftExons[len(leftExons)-1]:                                                    # check if lastLeftExon equals the exonNumber 
                            startNodeLastLeftExon = edgeKey[0]                                                                      # define startnode of last left exon
                        elif edgeValue['exon'] == rightExons[0]:                                                                # otherwise: check if firstRightExons equals the exonNumber 
                            endNodeFirstRightExon = edgeKey[1]                                                                      # define startnode of first right exon
                    if startNodeLastLeftExon!=-1 and endNodeFirstRightExon!=-1:                                                 # if startNode and endNode have been extracted skip the rest edges
                        break 
                full_path_dict = activeBinPathEnumeration3(str(endNodeFirstRightExon), str(startNodeLastLeftExon), [str(startNodeLastLeftExon)], {}, [0], [], graph, Bins) # employ pathEnumeration with MultiBinConstraint
                if full_path_dict != None:                                                                                      # If at least one path between real_last_left_exon and real_first_right_exon has been found
                    appendBinsLong(full_path_dict, bin, True, Bins, allNewBins, BinT, count, leftExons, rightExons)             # append it
    return (Bins+allNewBins), allNewBins                                                                                        # Return newMultiBins and allNewBins

    