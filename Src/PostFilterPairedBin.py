#!/usr/bin/env python

from importlib.resources import path
import sys, ast, os
from turtle import left
import networkx as nx
from collections import namedtuple
from sklearn.cluster import k_means

from sqlalchemy import false
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript
from PathEnumeration import activeBinPathEnumeration3
from PairedBinsToBins_short import checkPairedBins

def groupPairedBins(pairedBins):
    PairedBinT = namedtuple('GroupPairedBinT', 'leftExons rightExons count')
    groupedPairedBins = []
    predecessor = []
    rightExonList = []
    countList = []
    counter = 0
    for bin in pairedBins:
        # 1. Eliminate potential order mistakes
            # e. g. left   1,3,4    right: 5,6,7
        if checkPairedBins (bin) == False:
            print(bin)                                                                                              
            continue
        rightExons = bin.rightExons                                                                                                     # Assign bin.rightExons to rightExons
        leftExons = bin.leftExons                                                                                                       # Assign bin.leftExons to leftExons
        count = bin.count 
        
        # 2. Switch Cases to establish a general between right and left
        # e. g. left 1,2,5    right 0         ->  left 0          right 1,2,5
        if rightExons [0] < leftExons[0]:
            switchExons = rightExons
            rightExons = leftExons
            leftExons = switchExons 
        
        # 3. Eliminate cases, in which the first rightExon is not in leftExons, but smaller than last leftExons
        # e. g.   case 4    left 1,2,5      right 3,5,8
        if rightExons[0]<leftExons[len(leftExons)-1] and rightExons[0] not in leftExons:                                          
            continue

        # 4. Eliminate 
        # switched Case 11  left 0,1,2,8    right 1,2,5     -> left 0,1,2,8     right 5 -> discarded
        # and switched case 4.1:     right 1,2,5     left 2,3,5,8
        # No Triming
        j = 0
        incorrectPairedBinBoolean = False
        for i in range (0, len(leftExons)):                                                                                    
            if rightExons[j] == leftExons[i]:
                if j == len(rightExons)-1:
                    rightExons.pop(j-1)                                                                                        # Match: (e. g. case 11)
                    break
                j = j+1                                                                                                         # Remove firstRightExon
                if j < len(rightExons)-1:                                                                                       # If there's a remaining Exon in rightExons
                    if i<len(leftExons)-1: 
                        if leftExons[i+1] != rightExons[j]:                                                                     # check, whether there are additional elements in both, leftExons and rightExons                                                                                             # successors miss-match
                            incorrectPairedBinBoolean = True                                                                    # report it via incorrectpairedBinBoolean
                            break
                        elif leftExons[i+1] == rightExons[j]:
                            rightExons.pop(j-1)                                               
                else:                                                                                                           # If rightExons is empty exit for-loop
                    break                                                                                                       # Exit for-loop (since whileLoopBreakBoolean was initialized with False -> while-loop will not be exited)                                                                                                                                                    
            
        if incorrectPairedBinBoolean == True:                                                                                   # If an incorrect PairedBin has been found, discard it 
            continue
        
        # 5. Eliminate Bins with empty rightExons
        if len(rightExons) == 0:
            continue

        # 6. Eliminate Bins not resulting in real multBins, like
        # left: 0 right: 1
        # left: 2 right: 3
        # avoid GraphMistakes like going from 3->7->2 

        if len(leftExons)==1 and len(rightExons) == 1 and leftExons[0] + 1 == rightExons[0]:                                    # If we have case -1 (e. g. left 0, right 1)-> nothing to do
            continue     
        
        # 7. Group only valid and useful bins
        alreadyExistingBoolean = False
        if counter < 1:
            predecessor = leftExons
            rightExonList.append(rightExons)
            countList.append(count)
            counter = counter + 1
            continue
        if predecessor == leftExons:
            rightExonList.append(rightExons)
            countList.append(count)
        else:
            for bin in groupedPairedBins:
                if bin.leftExons==predecessor:
                    rightExonList = bin.rightExons + rightExonList
                    countList = bin.count + countList
                    groupedPairedBins.remove(bin)
                    groupedPairedBins.append(PairedBinT(leftExons=predecessor, rightExons=rightExonList, count=countList))    
                    alreadyExistingBoolean = True
            if alreadyExistingBoolean == False:
                groupedPairedBins.append(PairedBinT(leftExons=predecessor, rightExons=rightExonList, count=countList))
            
            predecessor = leftExons
            rightExonList = [rightExons]
            countList = [count]
            
    return groupedPairedBins

def eliminateInvalidPaths(paths:dict, groupedPairedBins, validPathNumber):
    for k in range(len(paths.keys())):
        path = paths[k]
        for bin in groupedPairedBins:                                                           # For every check every groupedBin
            compatibleLeftBinBoolean = False                                                    # Define two Booleans for finding of compatible leftBins 
            compatibleRightBinBoolean = False      
            for i in range(0, len(path)):                                                      # Iterate over all Exons in the path 
                if path[i] == bin.leftExons[0]:                                                # Position in the Path matches first Exon in leftExons
                    if len(bin.leftExons)>1:                                                    # If bin.leftExons is > 1, otherwise iteration over all Exons in the Bin is unnecessary
                        if len(path)-i>=len(bin.leftExons):                                    # Residual length of the current path is greater or equal than the length of the currently checked leftBin                                                                                        
                            for j in range(1, len(bin.leftExons)):                              # Now check every following Exon in bin.leftExons 
                                if path[i+j] != bin.leftExons[j]:                              # If you find a missmatch:
                                    break                                                       # break
                                elif j == len(bin.leftExons)-1:                                 # If you have reached the end, and found only matches 
                                    compatibleLeftBinBoolean = True                             # you foud a compatible leftBinBoolean
                            break                                                               # break from "for i in range (0, len(value)", anyway, because you won't find another starting spot
                        else:                                                                   # Residual length of the current path is smaller than len(bin.leftExons)
                            break                                                               # therefore -> not compatible leftBin expectable
                    else:                                                                       # Size of left.BinExon is 1 => Match 
                        compatibleLeftBinBoolean = True                                         # Report this match
                        break
                elif path[i] > bin.leftExons[0]:                                               # if current Exon of the path is greater than first Exon of bin.leftExons
                    break                                                                       # break also
            if compatibleLeftBinBoolean == True:                                                # Only if you have found a compatible leftBin
                for rightBin in bin.rightExons:                                                 # Iterate every Bin in bin.rightExons
                    for i in range(0, len(path)):                                              # Iterate every Exon in the Path 
                        if path[i] == rightBin[0]:                                             # If you have found a match = a potential starting point
                            if len(rightBin)>1:
                                if len(path)-i>=len(rightBin):                                 # Residual length of the current path is greater or equal than the length of currently checked rightBin
                                    for j in range(1, len(rightBin)):                           # Iterate every Exon of this particular rightBin
                                        if path[j+i] != rightBin[j]:                           # If you have found a mismatch
                                            break                                               # Break from "For j in range(1,len(rightBin)"
                                        elif j == len(rightBin)-1:                              # If you have reached the end of this rightBins
                                            
                                            compatibleRightBinBoolean = True                    # Report finding of a compatible Right Bin
                                    break                                                       # break from "for i in range...." anyway
                                else:                                                           # Residual length of the current path is smaller than the lenght of the currently checked rightBin
                                    break                                                       # therefore break from "for i in range...."
                            else:                                                               # Size of rightExon is 1 => Match
                                compatibleRightBinBoolean = True                                # report this match
                                break
                        elif path[i]>rightBin[0]:                                               # If currently checked exon in the current path is greater than the first Exon of this particular rightBin 
                            break                                                               # Break from "for i in range" -> stop this search
                    if compatibleRightBinBoolean == True:                                       # If a compatible rightBin has been found, stop searching
                        break                                                                   # and break from "for rightBin in bin.rightExons:"
            if compatibleLeftBinBoolean == True and compatibleRightBinBoolean == False:         # if both Booleans are correct -> 
                paths.pop(k)                                                                    # Remove current path
                validPathNumber[0] = validPathNumber[0] + 1
                break
    return paths                                                                                # return all filteredPaths
