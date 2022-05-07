#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript

# Full Path Enumeration

def fullPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, graph):  
    if v == '1':
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
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
        for j in range(1, len(bin[0])):                                                                                                                 # Suche nach Kompatiblen Bins
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

def activeBinPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, Bins:list):
    if v == '1':
        allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
        path_number[0] = path_number[0] + 1 
        return
    kantenTyp = None
    if len(pfad)>1:
        kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']
    if kantenTyp == 'SpliceJunction':
        activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']
        previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']
        compatibleBins, newActiveBins, SuccInMultiBinsList = getCompatibleBins(activeBins, activeExonNumber, previousExonNumber)                               
        if (len(compatibleBins)==0 and len(activeBins)>0):
            return
        if len(SuccInMultiBinsList) == 0:
            getNewActiveMultiBinsEasily(newActiveBins, Bins, activeExonNumber)
        else:
            getNewActiveBinsWithMulitBinRestriction(newActiveBins, Bins, activeExonNumber, SuccInMultiBinsList)        
    else:
        newActiveBins=activeBins
    for u in graph.adj[str(v)].keys():
        activeBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, Bins)
    return allpaths



def getMultiBins(Bins:list):
    multiBins = []
    for bin in Bins:
        if len(bin[0])>2:
            multiBins.append(bin)
    return multiBins

def getCompatibleMultiBins(act_Bins:list, actExon, prevExon):
    compatibleBins = []
    newActiveBins = []
    SuccInMultiBinsList = []
    for bin in act_Bins:
        for j in range(1, len(bin[0])):                                                                                                                 # Suche nach Kompatiblen Bins
            if bin[0][j] == actExon and bin[0][j-1] == prevExon:
                compatibleBins.append(bin)
                if j<len(bin[0])-1:
                    newActiveBins.append(bin)                                                                                                           # Nutze den Schleifendurchlaf, um aktive Multibins schon 
                                                                                                                                                        # in newActiveBins zu schreiben
                    SuccInMultiBinsList.append(bin[0][j+1])                                                                                             # Füge die Nachfolger des aktiven Exons in einem Multibin der 
                                                                                                                                                        # Liste SuccInMultiBinsList hinzu
    return compatibleBins, newActiveBins, SuccInMultiBinsList

def getNewActiveMultiBinsEasily (newActiveBins:list, Bins:list, act_Exon):
    for bin in Bins:    
        if (bin[0][0])==act_Exon:
            newActiveBins.append(bin)

def getNewActiveMultiBinsWithMultiBinRestriction(newActiveBins, Bins, act_Exon, SuccInMultiBinsList):
    for bin in Bins:
        if bin[0][0]==act_Exon and bin[0][1] in SuccInMultiBinsList:             
            newActiveBins.append(bin)

def activeMultiBinPathEnumeration(v:str, pfad:list, allpaths:dict, path_number:list, activeBins:list, graph, multiBins:list, counter:list):
    if len(multiBins)==0:
        return fullPathEnumeration(v, pfad, allpaths, path_number, graph)
    else:
        if v == '1':
            allpaths[path_number[0]] = nodepath_to_transcript(graph, pfad)
            path_number[0] = path_number[0] + 1 
            return
        kantenTyp = None
        if len(pfad)>1:
            kantenTyp = graph.edges[str(pfad[len(pfad)-2]), v]['type']
        if kantenTyp == 'SpliceJunction':
            activeExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['endExon']
            previousExonNumber=graph.edges[str(pfad[len(pfad)-2]), v]['startExon']
            compatibleBins, newActiveBins, SuccInMultiBinsList = getCompatibleMultiBins(activeBins, activeExonNumber, previousExonNumber)                               
            if len(compatibleBins)==0 and len(activeBins)>0:
                return
            if len(SuccInMultiBinsList) == 0:
                getNewActiveMultiBinsEasily(newActiveBins, multiBins, activeExonNumber)
            else:
                getNewActiveMultiBinsWithMultiBinRestriction(newActiveBins, multiBins, activeExonNumber, SuccInMultiBinsList)        
        else:
            newActiveBins=activeBins
        for u in graph.adj[str(v)].keys():
            activeMultiBinPathEnumeration(u, pfad+[u], allpaths, path_number, newActiveBins, graph, multiBins, counter)
    print(counter[0])
    return allpaths
    
