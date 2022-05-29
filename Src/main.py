import sys, ast, os
import networkx as nx
from collections import namedtuple
from PathEnumeration import fullPathEnumeration, activeBinPathEnumeration, activeBinPathEnumeration2, activeBinPathEnumeration3, getMultiBins
from PathEnumeration_old import activeMultiBinPathEnumeration
from PairedBinsToBins import fromPairedBinsToBins, checkPairedBins
from PostFilterPairedBin import groupPairedBins, eliminateInvalidPaths
from PairedBinsToBins_short import fromPairedBinsToBinsShort
from parse_graph_list_commented_Arbeitsdatei import parse_meta, parse_bins, parse_pairs, parse_graph, write_valid_gtf_entry, nodepath_to_transcript
from copy import deepcopy

dummyf = open("dummyout.gtf", "w")                                                                                              # output for the dummy code
dummyGeneCounter = 0

with open('human_geuvadis_simulated_5sets.graph') as f:
    fileEndReached = False
    f.readline()                                                                                                    #skip ---- seperator line
    
    # Define necessary Dictionaries                                                                                                     
    paths = {}
    filteredPaths = {}
    pairedBinDict = {}

    #Define necessary lists
    pfadListe = []

    # Define necessary Integer Varibales
    geneCounter = 0
    enumerationPathCounter = 0
    validPathCounter = 0
    
    while not fileEndReached:
        f.readline()                                                                                                # skip ==META: Read this line, but don't do anything
        Chromosome, Strand, Exons = parse_meta(f)                                                                   # Übergib f jetzt an def parse_meta, um Metadaten auszulesen und schreib diese in Chromosome, 
                                                                                                                    # Strand und Exons (Listen)
        Bins = parse_bins(f)                                                                                        # Lies die BINS aus f mit Hilfe der parse_bin(f) Funktion aus und schreib sie in Bins
        PairedBins = parse_pairs(f)
        
        G_full = nx.DiGraph()                                                                                       # Erzeug einen Diagraphen mit Hilfe von NetworkX
        fileEndReached, skip = parse_graph(f, G_full, Exons)                                                        # Ruf die Funktion parse_graph auf, übergibt ihr das File (f, den erzeugten Graphen und zugehörige Liste 
                                                                                                                    # mit den Exons), schreib den vollen Graphen in G_full, weise Skip einen Boolean zu, der true ist, wenn 
                                                                                                                    # die Zeile mit - beginnt                                                                                                                                                                                                      
        if not fileEndReached and not skip:                                                                         # Falls denoised Graph existiert, ruf wieder die Funktion Parse_Graph auf, übergib ihr das File 
                                                                                                                    # (f, Graph_clean und zugehörige Liste mit den Exons), schreib die # denoised Informationen aus dem 
                                                                                                                    # Graph-File in graph clean und übergib die letzte Zeile fileEndReached
            G_clean = nx.DiGraph()                                                                                  # Erzeug einen neuen gerichteten Graphen 
            fileEndReached, _ = parse_graph(f, G_clean, Exons)

        # Define necessary variables and assign
            # pathNumbers to enumerate paths  
        
        enumerationPathNumber = [0]
        validPathNumber = [0]
        pairedBinsCopy = deepcopy(PairedBins)

        # 1. Get MultiBins only from Bins
        MultiBins = getMultiBins(Bins)

        # 3. Add PairedBins to MultiBins
        MultiBins, gepaarteBins = fromPairedBinsToBinsShort(PairedBins, MultiBins, G_clean, Exons, True)
        
        # 4. Enumerate Paths with ActiveBinPathEnumeration3
        paths['Gene' + str(geneCounter)] = activeBinPathEnumeration3('1', '0', ['0'], {}, enumerationPathNumber, [], G_clean, MultiBins)

        # 5. Group paired Bins to filter valid paths after Enumeration 
        groupedPairedBins = groupPairedBins(pairedBinsCopy)

        # 6. Filter invalid Paths
        filteredPaths['Gene' + str(geneCounter)] = eliminateInvalidPaths(paths['Gene' + str(geneCounter)], groupedPairedBins, validPathNumber)
        
        # 7. Increase geneCounter
        geneCounter = geneCounter + 1

        # 8. Count paths
        enumerationPathCounter = enumerationPathCounter + enumerationPathNumber[0]
        validPathCounter = validPathCounter + validPathNumber[0]
        
        # All Paths Enumeration

        # Note: source and drain are ALWAYS named "0" and "1" respectively

        #if skip:
            # handle the rare case that noise deletion removes the whole second graph
        

        # TODO WORK WITH THE GRAPH HERE
        #Access Edge Types : G.edges[n1 , n2]['type'] == "Exon" || "SpliceJunction" || "Helper"
        #Access Main Coverage Count of an Edge : G.edges[n1 , n2]['counts']['c']
        #Access Exon length G.edges[n1 , n2]['length']
        
        #Source Node s is always G.nodes['0']
        #Drain Node t is always G.nodes['1']
        
        # DUMMY Code extracts longest Path (by number of bases) and writes it to a GTF file
        

        lpath = nx.dag_longest_path(G_full, weight="length")
        transcript = nodepath_to_transcript(G_full, lpath)
        write_valid_gtf_entry(dummyf, Chromosome, Strand, Exons, transcript, "Gene"+str(dummyGeneCounter), "Transcript"+str(dummyGeneCounter)+".1")
        dummyGeneCounter = dummyGeneCounter + 1

dummyf.close()
print(enumerationPathCounter)
print(validPathCounter)