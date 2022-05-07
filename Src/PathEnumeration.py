#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple

# Full Path Enumeration

def fullPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, graph):  
    if v == '1':
        allpaths[path_number[0]] = pfad
        path_number[0] = path_number[0] + 1
        return
    else:
        for u in graph.adj[str(v)].keys():
            fullPathEnumeration(u, pfad+[u], allpaths, path_number, graph)
    return allpaths

# Get Compatible Bins

def getCompatibleBins(act_Bins:list, actExon, prevExon):
    compatibleBins = []
    newActiveBins = []
    SuccInMultiBinsList = []
    for bin in act_Bins:
        for j in range(1, len(bin[0])-1):
            if bin[0][j] == actExon and bin[0][j-1] == prevExon:
                compatibleBins.append(bin)
                if j<len(bin[0])-1:
                    newActiveBins.append(bin)                                                                                                           # Nutze den Schleifendurchlaf, um aktive Multibins schon 
                                                                                                                                                        # in newActiveBins zu schreiben
                    SuccInMultiBinsList.append(bin[0][j+1])                                                                                             # Füge die Nachfolger des aktiven Exons in einem Multibin der 
                                                                                                                                                        # Liste SuccInMultiBinsList hinzu
    return compatibleBins, newActiveBins, SuccInMultiBinsList

# Get new Bins without MultiBinRestriction

def getNewActiveBinsEasily (newActiveBins:list, Bins:list, act_Exon):
    for bin in Bins:    
        if (bin[0][0])==act_Exon:
            newActiveBins.append(bin)

def getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, act_Exon, SuccInMultiBinsList):
    for bin in Bins:
        if bin[0][0]==act_Exon:
            if len(bin[0]) == 1: 
                newActiveBins.append(bin)
            else: 
                if bin[0][1] in SuccInMultiBinsList:             
                    newActiveBins.append(bin)

# Active Bin Enumeration

def activeBinPathEnumeration (v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v == '1':
        allpaths[path_number[0]] = pfad
        path_number[0] = path_number[0] + 1 
        return
    kantenTyp = None
    if len(pfad)>1:
        kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']
    if kantenTyp == 'SpliceJunction':
        activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']
        previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']
        compatibleBins, newActiveBins, SuccInMultiBinsList = getCompatibleBins(activeBins, activeExonNumber, previousExonNumber)                               
        if (len(compatibleBins)==0 and len(activeBins)==0):
            return
        if len(SuccInMultiBinsList) == 0:
            getNewActiveBinsEasily(newActiveBins, Bins, activeExonNumber)
        else:
            getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, activeExonNumber, SuccInMultiBinsList)        
    else:
        newActiveBins=activeBins
    for u in graph.adj[str(v)].keys():
        activeBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)
    return allpaths