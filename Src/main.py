# IMPORT
import sys
import parse_graph_new
import path_enumeration
import pairedbin_enumeration
import networkx as nx
import time
import matplotlib.pyplot as plt
from copy import deepcopy

start = time.time()

gene_id = 0
count = 1
no_trans = 0
#file_gtf = open("transcripts.gtf", "w") 
#textfile = "vergleich.txt"
#file2 = open(textfile, "w")

# MAIN
with open(sys.argv[1]) as f:
    fileEndReached = False
    f.readline() #skip ---- seperator line
    while not fileEndReached:
        
        # READ META AND BIN DATA FROM FILE
        f.readline() #skip ==META 
        Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
        Bins = parse_graph_new.parse_bins(f)
        PairedBins = parse_graph_new.parse_pairs(f)
        PairedBins_copy = deepcopy(PairedBins)

        # BUILD GRAPH
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons)

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()
            fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons)
            # nx.draw_networkx(G_clean, with_labels=True, arrowsize=12)
            # plt.show()
        if skip:
            G_clean = G_full

        transcripts = []
        
        # FULL PATH ENUMERATION
        #transcripts = path_enumeration.enumeration(G_clean,[],"0",["0"],"1")
        #no_trans = no_trans + len(transcripts)

        # BIN ENUMERATION
        #multi_bins = path_enumeration.get_multibins(Bins)
        #transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],multi_bins,"1")
        #no_trans = no_trans + len(transcripts)
        
        # PAIRED BIN ENUMERATION 1
        multi_bins = path_enumeration.get_multibins(Bins)
        paired_bins = pairedbin_enumeration.get_pairedbins(G_clean,PairedBins_copy,multi_bins)
        transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],paired_bins,"1")
        no_trans = no_trans + len(transcripts)
        
        # PAIRED BIN ENUMERATION 2
        #pairedbins_grouped = pairedbin_enumeration.group_pairs(PairedBins_copy)
        #transcripts = path_enumeration.enumeration_bins2(G_clean,[],"0",["0"],[],multi_bins,"1")
        #transcripts_copy = deepcopy(transcripts)
        #filtered_transcripts = pairedbin_enumeration.filter_transcripts(transcripts_copy,pairedbins_grouped)
        #no_trans = no_trans + len(filtered_transcripts)
        
        
        #laenge = len(filtered_transcripts)
        #file2.write("\n"+str(laenge))
        
        
        """
        if count == 695:
            
            print(PairedBins)
            print(pairedbins_grouped)
            print(transcripts)
            print(filtered_transcripts)
            
        count += 1
        """

        # ADD TRANSKRIPTS TO GTF FILE
        """
        gene_id += 1
        transcript_id = 0
        for transcript in transcripts:
            transcript_id += 1
            parse_graph_new.write_valid_gtf_entry(file_gtf,Chromosome,Strand,Exons,transcript,"Gene"+str(gene_id),"Transcript"+str(transcript_id))
        """
        
        #Access Edge Types : G.edges[n1 , n2]['type'] == "Exon" || "SpliceJunction" || "Helper"
        #Access Main Coverage Count of an Edge : G.edges[n1 , n2]['counts']['c']
        #Access Exon length G.edges[n1 , n2]['length']
        #Source Node s is always G.nodes['0']
        #Drain Node t is always G.nodes['1']
        # DUMMY Code extracts longest Path (by number of bases) and writes it to a GTF file
        
    end = time.time()
    
    print("Gesamtanzahl Transkripte: ", no_trans)
    # print("Anzahl paired bins:", no_bins)
    # print(liste)
    # print(len(liste))
    print('{:5.3f}s'.format(end-start))       
    #file_gtf.close()
    #file2.close()