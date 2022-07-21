#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple
import matplotlib as mpl
import matplotlib.pyplot as plt

ExonT = namedtuple('ExonT', 'leftPos rightPos')
BinT = namedtuple('BinT', 'exons count')
PairedBinT = namedtuple('PairedBinT', 'leftExons rightExons count')

# Data Description
Chromosome = None # String Name of the Chromosome this Gene is located on
Strand = None # String + or - as a descriptor of the strand of origin
Exons = None # List of ExonT
Bins = None # List of BinT
PairedBins = None # List of PairedBinT
G_full = None # The full splice-graph (as an NX Digraph/DAG)
G_clean = None # The splice-graph after basic cleaning steps (as an NX Digraph/DAG)

def parse_meta(f):
    Exons = []
    line = f.readline()
    Chromosome, Strand = line.strip().split()
    line = f.readline()
    while not line.startswith('='):
        splits = line.strip().split()
        Exons.append(ExonT(leftPos=int(splits[2]), rightPos=int(splits[3])))
        line = f.readline()

    return Chromosome, Strand, Exons

def parse_bins(f):
    Bins = []
    line = f.readline()
    while not line.startswith('='):
        splits = line.strip().split()

        exonIds = [ int(x) for x in splits[1].split(',') ]
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[3].split(';')[:-1] }
        
        Bins.append(BinT(exons=exonIds, count=counts))
        line = f.readline()
    return Bins

def parse_pairs(f):
    Bins = []
    line = f.readline()
    while not line.startswith('='):
        splits = line.strip().split()
        
        exonIdsLeft = [ int(x) for x in splits[1].split(',') ]
        exonIdsRight = [ int(x) for x in splits[3].split(',') ]        
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[5].split(';')[:-1] }
        
        Bins.append(PairedBinT(leftExons=exonIdsLeft, rightExons=exonIdsRight, count=counts))
        line = f.readline()
    return Bins

def parse_graph(f, G, Exons):
    line = f.readline()
    while not line.startswith("@arcs"):
        line = f.readline()
    f.readline() # skrip descriptor line

    line = f.readline() # first edgelist element
    while not line.startswith("@attributes"):
        splits = line.strip().split()
        
        leftId = splits[0]
        rightId = splits[1]
        edgetype, feature = splits[3].split(":")
        counts = { s.split(':')[0]:int(s.split(':')[1]) for s in splits[4].split(';')[:-1]   }
        
        if not G.has_node(leftId):
             G.add_node(leftId)
        if not G.has_node(rightId):
             G.add_node(rightId)
        
        G.add_edge(leftId , rightId)
        G.edges[leftId , rightId]['type'] = edgetype
        G.edges[leftId , rightId]['counts'] = counts
        G.edges[leftId , rightId]['length'] = 0 # default, overwritten for actual features
        if edgetype == "SpliceJunction":
            sExon, eExon = [int(i) for i in feature.split("-")]
            G.edges[leftId , rightId]['startExon'] = sExon
            G.edges[leftId , rightId]['endExon'] = eExon
            G.edges[leftId , rightId]['length'] = 1
        elif edgetype == "Exon":
            G.edges[leftId , rightId]['exon'] = int(feature)
            G.edges[leftId , rightId]['length'] = Exons[int(feature)].rightPos - Exons[int(feature)].leftPos + 1
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
    

# Workflow
"""
dummyf = open("dummyout.gtf", "w")  # output for the dummy code
dummyGeneCounter = 0

with open(sys.argv[1]) as f:
    fileEndReached = False
    f.readline() #skip ---- seperator line
    while not fileEndReached:
        f.readline() #skip ==META 
        Chromosome, Strand, Exons = parse_meta(f)
        print(Chromosome)
        Bins = parse_bins(f)
        print("Bins:", Bins)
        PairedBins = parse_pairs(f)
        print("Paired-Bins:", PairedBins)
        
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph(f, G_full, Exons)

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()
            fileEndReached, _ = parse_graph(f, G_clean, Exons)
            # nx.draw_networkx(G_clean, with_labels=True, arrowsize=12)
            # plt.show()
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
"""
