import sys, ast, os
import networkx as nx
from collections import namedtuple
from PathEnumeration import fullPathEnumeration, activeBinPathEnumeration, activeBinPathEnumeration2, activeBinPathEnumeration3, getMultiBins
from PathEnumeration_old import activeMultiBinPathEnumeration
from PairedBinsToBins import fromPairedBinsToBins, checkPairedBins
from PostFilterPairedBin import groupPairedBins
from PairedBinsToBins_short import fromPairedBinsToBinsShort
from parse_graph_list_commented_Arbeitsdatei import parse_meta, parse_bins, parse_pairs, parse_graph, write_valid_gtf_entry, nodepath_to_transcript

dummyf = open("dummyout.gtf", "w")                                                                                              # output for the dummy code
dummyGeneCounter = 0

with open('human_geuvadis_simulated_5sets.graph') as f:
    fileEndReached = False
    f.readline()                                                                                                     #skip ---- seperator line
    global_path_dict = {}
    #global_path_dict_multi = {}
    global_path_dict_full = {}
    gene_counter=0
    anzahlPfade = 0
    pairedBinDict = {}
    #anzahlPfade_multi=0                                                                                                   
    #anzahlPfade_full=0
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

        local_path_dict = {}
        #local_path_dict_multi={}
        #local_path_dict_full = {}
        pfad = ['0']
        #pfad_multi = ['0']
        pfad_full = ['0']
        path_number = [0]
        #path_number_multi = [0]
        #path_number_full = [0]
        zaehler = [0]
        MultiBins = getMultiBins(Bins)
        MultiBins=fromPairedBinsToBinsShort(PairedBins, MultiBins, G_clean, Exons, True)
        # groupedPairedBins = groupPairedBins(PairedBins)
        # MultiBins, newMultiBins = fromPairedBinsToBins (PairedBins, MultiBins, G_clean, Exons)
        # Bins = fromPairedBinsToBins (PairedBins, Bins, G_clean, Exons)
        # incorrectPairedBins = checkPairedBins(PairedBins)
        # Bins, newMultiBins = fromPairedBinsToBins (PairedBins, Bins, G_clean, Exons)
        #global_path_dict['Gene' + str(gene_counter)] = activeBinPathEnumeration3('0', pfad, local_path_dict, path_number, [], G_full, Bins)
        global_path_dict['Gene' + str(gene_counter)] = activeBinPathEnumeration3('1','0', pfad, local_path_dict, path_number, [], G_clean, MultiBins)        
        #global_path_dict_multi['Gene' + str(gene_counter)] = activeMultiBinPathEnumeration('0', pfad_multi, local_path_dict_multi, path_number_multi, MultiBins, G_clean, MultiBins)
        #global_path_dict_full['Gene' + str(gene_counter)] = fullPathEnumeration('0', pfad_full, local_path_dict_full, path_number_full, G_full)
        gene_counter = gene_counter + 1
        anzahlPfade = anzahlPfade + path_number[0]
        #anzahlPfade_full = anzahlPfade_full + path_number_full[0]
        #anzahlPfade_multi = anzahlPfade_multi + path_number_multi[0]
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
print(anzahlPfade)