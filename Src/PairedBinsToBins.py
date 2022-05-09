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
    counter = 0
    for bin in pairedBins: 
        if counter == 0:
            addedBins = []
            counter = 1
        real_last_left_exon = bin.leftExons[len(bin.leftExons)-1]
        real_first_right_exon = bin.rightExons[0]

        last_left_exon = bin.leftExons[len(bin.leftExons)-1]
        first_right_exon = bin.rightExons[0]
        while(last_left_exon >= first_right_exon and (len(bin.leftExons)>1 or len(bin.rightExons)>1)):
            if len(bin.leftExons)>1:
                bin.leftExons.remove(bin.leftExons[len(bin.leftExons)-1])
                last_left_exon = bin.leftExons[len(bin.leftExons)-1]
            else:
                if len(bin.rightExons)>1:
                    bin.rightExons.remove(bin.rightExons[0])
                    first_right_exon = bin.rightExons[0]
        exonBin = bin.leftExons + bin.rightExons
        

        # if (len(bin.leftExons) + len(bin.rightExons)==2 and int(last_left_exon) + 1 == int(first_right_exon)):
        #     break 
        # if last_left_exon == first_right_exon:
        #     if len(bin.leftExons)>=len(bin.rightExons):
        #         exons = bin.leftExons + bin.rightExons[0:len(bin.rightExons)-2]

        # if first_right_exon >= last_left_exon:
        if exonBin not in addedBins and real_last_left_exon>=real_first_right_exon:
            startNodeLastLeftExon = real_last_left_exon + 2
            endNodeFirstRightExon = len(Exons)*2+1-real_first_right_exon 
            full_path_dict = PairedBinToBinEnumeration(str(startNodeLastLeftExon), str(endNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], graph)
            if full_path_dict != None:
                for key, value in full_path_dict.items():
                    if len(value)>2:
                        newBin = bin.leftExons + value[1:len(value)-2] + bin.rightExons
                    else:
                        newBin = bin.leftExons + bin.rightExons
                    newBinBoolean = True
                    for bin1 in Bins:
                        if newBin == bin1.exons:
                            newBinBoolean = False
                            if counter > 1:
                                break
                             
                    # Add all new Bins implied from Paired Bins to Bins
                    addedBins.append(newBin)
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
