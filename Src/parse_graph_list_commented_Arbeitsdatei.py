#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple

ExonT = namedtuple('ExonT', 'leftPos rightPos')
BinT = namedtuple('BinT', 'exons count')
PairedBinT = namedtuple('PairedBinT', 'leftExons rightExons count')

# Data Description
Chromosome = None                                                                                                   # String Name of the Chromosome this Gene is located on
Strand = None                                                                                                       # String + or - as a descriptor of the strand of origin
Exons = None                                                                                                        # List of ExonT
Bins = None                                                                                                         # List of BinT
PairedBins = None                                                                                                   # List of PairedBinT
G_full = None                                                                                                       # The full splice-graph (as an NX Digraph/DAG)
G_clean = None                                                                                                      # The splice-graph after basic cleaning steps (as an NX Digraph/DAG)

def parse_meta(f):                                                                                                  # Read line aus file f (damit sind wir jetzt in Zeile 3, in der "GL00008.2 steht)
    Exons = []                                                                                                      # Schreib 1. Teil des Strings (GL00008.2) bis zum Leerzeichen in Chromosome und 2. Teil (-) in die Variable 
                                                                                                                    # Strand, in # dem erst alle white spaces entfernt werden und dann der String am Leerzeichen getrennt wird
    line = f.readline()                                                                                             # Lies nächste Zeile des Graphfiles (z. B. "Exon 0 179629 179806")
    Chromosome, Strand = line.strip().split()                                                   
    line = f.readline()
    while not line.startswith('='):                                                                                 # Bis ==Bins kommt tue
        splits = line.strip().split()                                                                               # Teil den String in der Zeile anhand der Leerzeichen und schreib jeden String einzeln in eine Liste mit dem 
                                                                                                                    # Namen Split 
        Exons.append(ExonT(leftPos=int(splits[2]), rightPos=int(splits[3])))                                        # Füg der Liste an Exons das Tupel aus ExonT, leftPos und rightPos hinzu
        line = f.readline()                                                                                         # Wiederhol das ganze für alle Zeilen, in denen die Exons aufgelistet sind

    return Chromosome, Strand, Exons                                                            

# BINS repräsentieren die möglichen Verknüfpungen der Exons miteinander, 
# die aus den Informationen der jeweiligen einzelnen Reads gezogen werden

def parse_bins(f):                                                                                                  # Funktion zum Auslesen der Daten der BINS 
    Bins = []                                                                                                       # Erzeug leeres Array/Liste für die Bins
    line = f.readline()                                                                                             # Lies neue Zeile, sodass wir jetzt bei "Bin 0,1 Count ...." sind
    while not line.startswith('='):                                                                                 # Solange wir nicht bei den Paired Bins angelangt sind
        splits = line.strip().split()                                                                               # Teil den String "line" anhand der Leerzeichen und schreib ihn in die Liste splits

        exonIds = [ int(x) for x in splits[1].split(',') ]                                                          # Weis der Variable ExonID alle Werte (getrennt nach Komma: z. B. 1,2 => 1 und 2) im zweiten Element 
                                                                                                                    # der Split-Liste zu
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[3].split(';')[:-1] }                        # Trenn das 4. Split-Listenelement nach ";" und ignorier das letzte (leere) Element dieser Liste. 
                                                                                                                    # Trenn jedes Element dieser entstandenden Liste wiederum nach ":" und füg dem Dictionary "counts" die 
                                                                                                                    # counts so hinzu, dass der 1. Teil der - die Nummer des Experiments als String (z. B. '1' bzw. am Ende 'c') 
                                                                                                                    # - der Key ist und der 2. Teil - die Anzahl an Counts in diesem # Experiment - der Value als Integer für 
                                                                                                                    # diesen Bin
                                                                                                                    
                                                                                                
                                                                                                
        Bins.append(BinT(exons=exonIds, count=counts))                                                              # Füg der Liste "Bins" das Tupel bestehend aus der Split-Junction (als Liste z. B. [0,1]) und den Counts 
                                                                                                                    # als Dictionary (z. B.
        line = f.readline()                                                                                         # Lies nächste Zeile solange bis wir bei den Paired Bins sind
    return Bins

# Paired BINS repräsentieren die möglichen Verknüpfungen der Exons miteinander, die aus den
# Informationen der jeweiligen Paired-Reads gezogen werden (z. B. 0-1 links und 2 rechts)

def parse_pairs(f):                                                                                                 # Funktion zum Auslesen der Daten über die Paired BINS       
    Bins = []                                                                                                       # Creat new empty Array for Bins
    line = f.readline()                                                                                             # Lies die nächste Zeile um in der 1. Zeile mit den Infomraitonen über die Paired Bins zu landen
    while not line.startswith('='):                                                                                 # Solange die Zeile nicht mit = beginnt (d. h. wir im Graphenteil des Skripts angekommen sind)
        splits = line.strip().split()                                                                               # Teil den String "line" anhand der Leerzeichen und schreib ihn in die Liste splits
        
        exonIdsLeft = [ int(x) for x in splits[1].split(',') ]                                                      # Schreib die Ids aller Exons (getrennt nach ,) vom linken Read (read1 in FW Richtung) aus dem 2.  
                                                                                                                    # Listenelement in eine Liste
        exonIdsRight = [ int(x) for x in splits[3].split(',') ]                                                     # Schreib die Ids aller Exons (getrennt nach ,) vom rechten Read (read2 in Rev-Richtung) aus dem 
                                                                                                                    # 4. Listenelement in eine Liste
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[5].split(';')[:-1] }                        # Trenne das 6. Split-Listenelement nach ";" und ignoriere das letzte (leere) Element dieser Liste. 
                                                                                                                    # Trenne jedes Element dieser # entstandenden Liste wiederum nach ":" und füg dem Dictionary "counts" 
                                                                                                                    # die counts so hinzu, dass der 1. Teil - die Nummer des Experiments als String (z. B. '1' bzw. am Ende 'c')
                                                                                                                    #  - der Key ist und der 2. Teil - die Anzahl an Counts in diesem Experiment - der Value als Integer für 
                                                                                                                    # diesen Bin
                                                                                                                    
        
        Bins.append(PairedBinT(leftExons=exonIdsLeft, rightExons=exonIdsRight, count=counts))                       # Füg der Liste Bins das Tupel bestehend aus den linken und rechten Split-Junctions (als Liste z. B. [0,1] [2]) und den Counts als 
                                                                                                                    # Dictionary (z. B. Experiment "1":1 und Gesamtcount"c":4) hinzu)
        line = f.readline()                                                                                         # Lies nächste Zeile solange bis wir bei den Paired Bins sind
    return Bins

def parse_graph(f, G, Exons):                                                                   
    line = f.readline()                                                                                             # Skip die 1. Zeile ("== Graph")
    while not line.startswith("@arcs"):                                                                             # Skip all lines until @arcs (including @arcs) and don't do anything
        line = f.readline()                                                                     
    f.readline()                                                                                                    # Skip descriptor line    

    line = f.readline()                                                                                             # first edgelist element
    while not line.startswith("@attributes"):                                                                       # Until the end of the table is not reached, do:
        splits = line.strip().split()                                                                               # Split lines according to spaces and write them into the list "splits"
        
    # Erklärungen zum Aufbau des Graphen den Id's:
        # Grundsätzlich: 
            # Exons werden nicht als Knoten sondern als Kanten betrachtet. 
            # Quelle und Target sind aber Knoten
            # ID (Quelle) := 0
            # ID (Target) := 1
            # Helperkanten sind Kanten aus der Quelle und in die Senke

            # Für Quellenkanten gilt:
                # leftID(Quellekante) := 0 
                # rightID(Quellenkante) := leftID(ExonZuDemSieFührt)
            # Für Senkenkanten gilt:
                # leftID(Senken) := rightID(ExonVonDemSieHerkommt)
                # rightID(Senkenkante) : = 2 
            # Für Kanten der Exons gilt:
                # leftID(Exon_i) : = i+2 . 
                # rightID(Exon_i) := #Exons*2+1-i

        leftId = splits[0]                                                                                          # Lies die linke Id aus dem 1. Listenelement von Splits 
        rightId = splits[1]                                                                                         # Lies die rechte Id aus dem 2. Listenelement von Splits
        edgetype, feature = splits[3].split(":")                                                                    # Aus dem 4. Listenelement lies den Typ des Knotens aus (z.B. SpliceJunction, Helper oder Exon) und 
                                                                                                                    # schreib ihn in die Varibale edgetype, lies as "feature" an sich aus (z. B. 0-1 für SpliceJunction,
                                                                                                                    #  None für Helper und 1 für Exon) und weis der Variable feature zu
                                                                                                
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[4].split(';')[:-1]   }                      # Trenn das 5. Split-Listenelement nach ";" und ignoriere das letzte (leere) Element dieser Liste. 
                                                                                                                    # Trenn jedes Element dieser # entstandenden Liste wiederum nach ":" und füg dem Dictionary "counts" 
                                                                                                                    # die counts so hinzu, dass der 1. Teil - die Nummer des # Experiments als String (z. B. '1' bzw. am Ende 'c') 
                                                                                                                    # - der Key ist und der 2. Teil - die Anzahl an Counts in diesem Experiment - # der Value als Integer für 
                                                                                                                    # diesen Graphen
                                                                                                                    
        
        # Add the nodes, no doublings

        if not G.has_node(leftId):                                                              
             G.add_node(leftId)
        if not G.has_node(rightId):
             G.add_node(rightId)
        
        # Add edges (ID rules: see above)

        G.add_edge(leftId , rightId)                                                                                # Add edge 
        G.edges[leftId , rightId]['type'] = edgetype                                                                # Write the feature edgetype (Splice-junction, Exon, Helper) into type
        G.edges[leftId , rightId]['counts'] = counts                                                                # Add the dictionary with the counts
        G.edges[leftId , rightId]['length'] = 0                                                                     # default, overwritten for actual features
        if edgetype == "SpliceJunction":
            sExon, eExon = [int(i) for i in feature.split("-")]                                                     # If it´s a splice junction
            G.edges[leftId , rightId]['startExon'] = sExon                                                          # Add First Exon 
            G.edges[leftId , rightId]['endExon'] = eExon                                                            # Add second Exon   
            G.edges[leftId , rightId]['length'] = 1                                                                 # length is 1
        elif edgetype == "Exon":                                                                                    # For exons
            G.edges[leftId , rightId]['exon'] = int(feature)                                                        # add exon number
            G.edges[leftId , rightId]['length'] = Exons[int(feature)].rightPos - Exons[int(feature)].leftPos + 1    # Add Exon length
        # else  edgetype == "Helper":
        
        line = f.readline()

    
    f.readline()
    f.readline()
    line = f.readline()
    
    return line == "", line.startswith('-')


# f ..  file to write to
# chrom  ..  string of the chromosome as produced by the parser
# strand  ..  string of the strand as produced by the parser
# exons  ..  all exons of this gene as list of ExonT as produced by the parser
# transcript  .. ordered list of indices of visited exons e.g. [0, 2, 3, 5]
# gene_id  ..  identifier of the current locus potentially containing multiple isoforms
# transcript_id  ..  identifier of the transcript
def write_valid_gtf_entry(f, chrom, strand, exons, transcript, geneId, transcriptId):

    # write transcript group line
    f.write(chrom+"\tFortgMethoden\ttranscript\t"+str(exons[transcript[0]].leftPos)+"\t"+str(exons[transcript[-1]].rightPos)+"\t0\t"+strand+'\t.\tgene_id "'+geneId+'"; transcript_id "'+transcriptId+'";\n')

    # write each exon in order, directly adjacent exons need to be joined
    
    left = exons[0].leftPos
    right = exons[0].rightPos
    for ti in transcript[1:]:
        if exons[ti].leftPos == right + 1: # adjacent, extend
            right = exons[ti].rightPos
        else:
            f.write(chrom+"\tFortgMethoden\texon\t"+str(left)+"\t"+str(right)+"\t0\t"+strand+'\t.\tgene_id "'+geneId+'"; transcript_id "'+transcriptId+'";\n')    
            left = exons[ti].leftPos
            right = exons[ti].rightPos 
                   
    f.write(chrom+"\tFortgMethoden\texon\t"+str(left)+"\t"+str(right)+"\t0\t"+strand+'\t.\tgene_id "'+geneId+'"; transcript_id "'+transcriptId+'";\n')    

def nodepath_to_transcript(G, path):
    
    transcript = []
    for n1, n2 in zip(path, path[1:]):
        if G.edges[n1 , n2]['type'] == "Exon":
            transcript.append( G.edges[n1 , n2]['exon'] )
    return transcript



""" dummyf = open("dummyout.gtf", "w")                                                                                              # output for the dummy code
dummyGeneCounter = 0

with open('test.graph') as f:
    fileEndReached = False
    f.readline()                                                                                                     #skip ---- seperator line
    global_path_dict = {}
    global_path_dict_full = {}
    gene_counter=0
    anzahlPfade = 0                                                                                                   
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
        local_path_dict_full = {}
        pfad = ['0']
        pfad_full = ['0']
        path_number = [0]
        path_number_full = [0]   

        # Bins = fromPairedBinsToBins (PairedBins, Bins, G_clean)
        global_path_dict['Gene' + str(gene_counter)] = activeBinPathEnumeration('0', pfad, local_path_dict, path_number, Bins, G_full, Bins)
        
        # global_path_dict_full['Gene' + str(gene_counter)] = Enumeration.fullPathEnumeration('0', pfad_full, local_path_dict_full, path_number_full)
        #gene_counter = gene_counter + 1
        anzahlPfade = anzahlPfade + path_number[0]
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
print(global_path_dict)
 """