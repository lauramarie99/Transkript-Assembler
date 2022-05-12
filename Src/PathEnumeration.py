#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple

from sqlalchemy import false
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
def getNewActiveBinsEasily (newActiveBins:list, Bins:list, current_Exon):
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

#Enumerate Paths with for-Loop later
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
        secondBinListBoolean = False
        if len(activeBins)!=0:
            for bin in activeBins:
                for j in range(1, len(bin.exons)):                                                                                                                 # Suche nach Kompatiblen Bins
                    if bin.exons[j] == activeExonNumber and bin.exons[j-1] == previousExonNumber:
                        # There#s at least one remaing active Bin
                        activeBinListBoolean = True
                        if j<len(bin.exons)-2:
                            newActiveBins.append(bin)                                                                                                           # Nutze den Schleifendurchlaf, um aktive Multibins schon                                                                                     # in newActiveBins zu schreiben
                            SuccInMultiBinsList.append(bin.exons[j-1])                                                                                             # Füge die Nachfolger des aktiven Exons in einem Multibin der 
                        elif bin.exons[j] == previousExonNumber:
                            secondBindListBoolean = True
                                                                                                                                            # Liste SuccInMultiBinsList hinzu
            if activeBinListBoolean == False and secondBinListBoolean == True:
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

#Enumerate Paths with ForLoop initially
def activeBinPathEnumeration2(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v == '1':
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
        path_number[0] = path_number[0] + 1 
        return
    newActiveBins = []
    for u in graph.adj[str(v)]:
        if graph.edges[v, u]['type'] == 'SpliceJunction':
            activeExonNumber=graph.edges[v, u]['endExon']
            previousExonNumber=graph.edges[v,u]['startExon']
            SuccInMultiBinsList = []
            activeBinListBoolean = False
            secondBinListBoolean = False
            for bin in activeBins:
                for j in range(0, len(bin.exons)-1):
                    if bin.exons[j+1] == activeExonNumber and bin.exons[j] == previousExonNumber:    
                        if j!=0:
                            activeBinListBoolean = True
                        if j!=len(bin.exons)-2:
                            newActiveBins.append(bin)
                    elif (bin.exons[j] == previousExonNumber) and j!=0:
                        secondBinListBoolean = True
                        SuccInMultiBinsList.append(bin.exons[j])                                                            # Füge die Nachfolger des aktiven Exons in einem Multibin der
                                                                                                           
            if activeBinListBoolean == False and secondBinListBoolean == True:
                return
            for bin1 in Bins:
                if bin1.exons[0] == previousExonNumber and bin1.exons[1] == activeExonNumber:
                    newActiveBins.append(bin1)
            #if len(SuccInMultiBinsList) == 0:
            #    getNewActiveBinsEasily(newActiveBins, Bins, activeExonNumber)
            #else:
            #    getNewActiveMultiBinsWithMultiBinRestriction(newActiveBins, Bins, activeExonNumber, SuccInMultiBinsList)        
        else:
            newActiveBins=activeBins
        activeBinPathEnumeration2(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)
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
        if bin.exons[0]==act_Exon and len(bin.exons)>2 and bin.exons[1] in SuccInMultiBinsList:             
            newActiveBins.append(bin)