#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript

# Full Path Enumeration

def fullPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, graph):  
    # If current is drain, return current path
    if v == '1':
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
        path_number[0] = path_number[0] + 1
        return
    else:
        #For all successors of the current node repeat this procedure
        for u in graph.adj[str(v)].keys():
            fullPathEnumeration(u, pfad+[u], allpaths, path_number, graph)
    return allpaths

   
# Get new Bins without MultiBinRestriction
def getNewActiveBinsEasily (newActiveBins:list, Bins:list, act_Exon):
    for bin in Bins:    
        if (bin.exons[0])==act_Exon:
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

def activeBinPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    # If current node ist drain, return current path
    if v == '1':
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
        path_number[0] = path_number[0] + 1 
        return
    kantenTyp = None
    # If pathlength is great than 1, retrieve ed
    if len(pfad)>1:
        kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']
    # Only if edgeType is a SpliceJunction
    if kantenTyp == 'SpliceJunction':
        activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']
        previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']
        newActiveBins = []
        SuccInMultiBinsList = []
        activeBinListBoolean = False
        if len(activeBins)!=0:
            for bin in activeBins:
                for j in range(1, len(bin.exons)):                                                                                                                 # Suche nach Kompatiblen Bins
                    if bin.exons[j] == activeExonNumber and bin.exons[j-1] == previousExonNumber:
                        # There#s at least one remaing active Bin
                        activeBinListBoolean = True
                        if j<len(bin.exons)-1:
                            newActiveBins.append(bin)                                                                                                           # Nutze den Schleifendurchlaf, um aktive Multibins schon 
                                                                                                                                                        # in newActiveBins zu schreiben
                            SuccInMultiBinsList.append(bin.exons[j+1])                                                                                             # Füge die Nachfolger des aktiven Exons in einem Multibin der 
                                                                                                                                                        # Liste SuccInMultiBinsList hinzu
            if activeBinListBoolean == False:
                return
        # Skip path, if activeBinList was not empty upon entering this node, but is now 
        if len(SuccInMultiBinsList) == 0:
            getNewActiveBinsEasily(newActiveBins, Bins, activeExonNumber)
        else:
            getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, activeExonNumber, SuccInMultiBinsList)        
    # If edgeType is not a SpliceJunction
    else:
        newActiveBins=activeBins
    for u in graph.adj[str(v)].keys():
        activeBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)
    return allpaths


#Get only Bins with length>2
def getMultiBins(Bins:list):
    multiBins = []
    for bin in Bins:
        if len(bin.exons)>2:
            multiBins.append(bin)
    return multiBins

def getNewActiveMultiBinsWithMultiBinRestriction(newActiveBins, multiBins, act_Exon, SuccInMultiBinsList):
    for bin in multiBins:
        if bin.exons[0]==act_Exon and len(bin.exons)>1 and bin.exons[1] in SuccInMultiBinsList:             
            newActiveBins.append(bin)

#Enumerate Paths with MultBins only 
def activeMultiBinPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, multiBins:list):
    # If the set containing multibins only is empty use full-Pathenumeration
    if len(multiBins)==0:
        return fullPathEnumeration(v, pfad, allpaths, path_number, graph)
    # If the set is not empty
    else:
        # If current Node is drain, return current path
        if v == '1':
            allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
            path_number[0] = path_number[0] + 1 
            return
        kantenTyp = None
        # If pathlength is greater then 1, define edgeType 
        if len(pfad)>1:
            kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']
        # If EdgeType is Splicejunction
        if kantenTyp == 'SpliceJunction':
            activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']
            previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']
            # Check for compatibleBins
            newActiveBins = []
            SuccInMultiBinsList = []
            activeBinListBoolean = False
            if len(activeBins)!=0:
                # Check every Bin in currentsBins
                for bin in activeBins:
                    for j in range(1, len(bin.exons)):                                                                                                              # Suche nach Kompatiblen Bins
                        # If current edge is represented in this Bin add bin to compatible Bins
                        if bin.exons[j] == activeExonNumber and bin.exons[j-1] == previousExonNumber:
                            activeBinListBoolean = True                                                                                                             # There's at least one remaining active Bin
                            # If the edge is not represented as the last two Exons in this Bin
                            if j<len(bin.exons)-1:
                                #Add Bin to NewActive Bins
                                newActiveBins.append(bin)                                                                                                           # Nutze den Schleifendurchlaf, um aktive Multibins schon 
                                #Add the Exon following the current edge in this Bin to SuccInMultiBinList                                                                                                                                      # in newActiveBins zu schreiben
                                SuccInMultiBinsList.append(bin.exons[j+1])                                                                                          # Füge die Nachfolger des aktiven Exons in einem Multibin der 
                # Skip path, if activeBinList was not empty upon entering this node, but is now                                                                     # Liste SuccInMultiBinsList hinzu                              
                if activeBinListBoolean == False:
                    return
            # If there's no Exon following this edge represented by an active MultiBin
            if len(SuccInMultiBinsList) == 0:
                getNewActiveBinsEasily(newActiveBins, multiBins, activeExonNumber)
            # If there's an Exon following this edge represented by an active MultiBin
            else:
                getNewActiveMultiBinsWithMultiBinRestriction(newActiveBins, multiBins, activeExonNumber, SuccInMultiBinsList)        
        # If edgeType is not a SpliceJunction
        else:
            newActiveBins=activeBins
        for u in graph.adj[str(v)].keys():
            activeMultiBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, multiBins)
    return allpaths
    
