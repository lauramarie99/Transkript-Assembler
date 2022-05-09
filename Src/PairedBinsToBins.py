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

def fromPairedBinsToBins(pairedBins, Bins, graph, Exons:list):
    BinT = namedtuple('BinT', 'exons count')
    # Lies Bins aus PairedBins
    for bin in pairedBins: 
        last_left_exon = bin.leftExons[len(bin.leftExons)-1]
        first_right_exon = bin.rightExons[0]

        if (len(bin.leftExons)+ len(bin.rightExons)==2 and int(last_left_exon) + 1 == int(first_right_exon)):
            break
        if first_right_exon >= last_left_exon:
            
            startNodeLastLeftExon = last_left_exon + 2
            endNodeFirstRightExon = len(Exons)*2+1-first_right_exon 
            full_path_dict = PairedBinToBinEnumeration(str(startNodeLastLeftExon), str(endNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], graph)
            if full_path_dict != None:
                for key, value in full_path_dict.items():
                    newBin = bin.leftExons + value + bin.rightExons
                    newBinBoolean = True
                    for bin1 in Bins:
                        if newBin == bin1.exons:
                            newBinBoolean = False
                            break 
                    if newBinBoolean==True:
                        Bins.append(BinT(exons=newBin, count={}))
    return Bins
            
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
