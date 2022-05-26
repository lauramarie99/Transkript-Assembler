import Parsing
import networkx as nx
from matplotlib import pyplot as plt
def get_bins(bins: list, new_bins: list, exon: int):
    for bin in bins:
        if bin.exons[0] == exon:
            new_bins.append(bin)
    return new_bins

successor = []
def MyEnu(Graph, Transcript, Path, Node):
    #Node is a drain if its "1"
    if Node == "1":
        Transcript.append(Parsing.nodepath_to_transcript(Graph,Path))
        return
    else:
        #Graph adj object holding neighbors of each node
        successor = Graph.adj[Node]
        for node in successor:
            MyEnu(Graph,Transcript, (Path + [node]), node)

        return Transcript
        print("result:")
        print(Transcript)


#working progress, does not work correctly yet
def EnuBinConstraint (Graph, Transcript, Path, Node, Bin):
    if Node == "1":
        Transcript.append(Parsing.nodepath_to_transcript(Graph,Path))
        return
    else:
        #Graph adj object holding neighbors of each node
        successor = Graph.adj[Node]
        for node in successor:
            ActiveBin = Bin
            if Graph.edges[Node, node]['type'] == "SpliceJunction":
                StartExon = Graph.edges[Node, node]["startExon"]
                EndExon = Graph.edges[Node, node]["endExon"]
                #get compatible bins
                for bin in ActiveBin:
                    for i in range(0, (len(bin.exons) - 1)):
                        if(StartExon != bin.exons[i] and EndExon != bin.exons[i+1]):
                            ActiveBin.remove(bin)
            #if ActiveBins now empty stop it
            if(len(ActiveBin) == 0):
                continue
            else:
                 EnuBinConstraint(Graph,Transcript, (Path + [node]), node, Bin)

    return Transcript 


#code that acutally works (from laura) 
def EnuBinConstraint(Graph, Transcript, Path, Node, bins, Bins,NodeNr):
    if Node == NodeNr:
        Transcript.append(Parsing.nodepath_to_transcript(Graph, Path))
        return
    else:
        successor = Graph.adj[Node]
        for node in successor:
            if Graph.edges[Node, node]['type'] == "SpliceJunction":
                # print("I AM HERE")
                StartExon = Graph.edges[Node, node]["startExon"]
                EndExon = Graph.edges[Node, node]["endExon"]
                compatible = False
                copy = []
                for bin in bins:
                    print(bin)
                    for i in range(0, (len(bin.exons) - 1)):
                        if ((bin.exons[i] == StartExon) and (bin.exons[i + 1] == EndExon)):
                            compatible = True
                            if (i != (len(bin.exons) - 2)):
                                copy.append(bin)

                if (compatible == False and (len(bins) != 0)):
                    continue

                for bim in Bins:
                        if ((bim.exons[0] == StartExon) and (bim.exons[1] == EndExon)):
                                copy.append(bim)
            else:
                copy = bins

            EnuBinConstraint(Graph, Transcript, Path + [node], node, copy, Bins,NodeNr)

    return Transcript



