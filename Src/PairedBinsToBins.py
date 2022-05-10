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
    counter = 0
    # Read Bins from PairedBins
    for bin in pairedBins: 
        #If this the first ListElement of PairedBins, create a List, to save all Bins, that have already been added
        if counter == 0:
            addedBins = []
            counter = 1
        # Define LastExon of left Bins and FirstExon of RightBins
        real_last_left_exon = bin.leftExons[len(bin.leftExons)-1]
        real_first_right_exon = bin.rightExons[0]
        # Make a working-copy of these Exons 
        last_left_exon = real_last_left_exon
        first_right_exon = real_first_right_exon
        #As long as last_left_exon is greater or equal to first_right_exon and there are at least two Elements in leftExons or rightExons
        while(last_left_exon >= first_right_exon and (len(bin.leftExons)>1 or len(bin.rightExons)>1)):
            #If there at least two exons in leftExons
            if len(bin.leftExons)>1:
                #Remove the last exon of leftExons and assign the new last exon of leftExons to the variable last_left_exon
                bin.leftExons.remove(bin.leftExons[len(bin.leftExons)-1])
                last_left_exon = bin.leftExons[len(bin.leftExons)-1]
            # If there's only one Exon left in leftExons
            else:
                # If there's more than one exon left in rightExons
                if len(bin.rightExons)>1:
                    #Remove the first eExon of rightExons and assign the new first exon of rightExons to the variable first_right_exon
                    bin.rightExons.remove(bin.rightExons[0])
                    first_right_exon = bin.rightExons[0]
        # Write the remaing Exons in leftExons and rightExons in exonBin
        exonBin = bin.leftExons + bin.rightExons

        # If this particular exonBin has not been added to addedBin and real_last_left_exon>=realfirst_right_exon
        if exonBin not in addedBins and real_last_left_exon>=real_first_right_exon:
            startNodeLastLeftExon = real_last_left_exon + 2
            endNodeFirstRightExon = len(Exons)*2+1-real_first_right_exon 
            full_path_dict = PairedBinToBinEnumeration(str(startNodeLastLeftExon), str(endNodeFirstRightExon), [str(startNodeLastLeftExon)], {}, [0], graph)
            #If at least one path between real_last_left_exon and real_first_right_exon has been found
            if full_path_dict != None:
                for value in full_path_dict.values():
                    if len(value)>2:
                        newBin = bin.leftExons + value[1:len(value)-2] + bin.rightExons
                    else:
                        newBin = bin.leftExons + bin.rightExons
                    newBinBoolean = True
                    for bin1 in Bins:
                        if newBin == bin1.exons:
                            newBinBoolean = False
                            break
                    # Add all new Bins implied from PairedBins to AddedBins
                    addedBins.append(newBin)
                    # If this Bin is not already in Bins
                    if newBin and newBinBoolean==True:
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
