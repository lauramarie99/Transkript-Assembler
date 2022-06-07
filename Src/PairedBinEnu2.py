# quick fix for personal script that uses elements from laura
import FinalEnu as path_enumeration
from collections import namedtuple

BinT = namedtuple('BinT', 'exons count')
def get_pairedbins(graph, pairedbins, multibins):
    all_pairedbins = []  # List storing all bins
    for pairedbin in pairedbins:
        start_node = ""  # Start node of enumeration
        end_node = ""  # End node of enumeration
        new_pairedbins = []  # New bins, representing all possible paths from the left to the right exons of each paired bin
        left = pairedbin.leftExons  # Left exons of paired bin
        right = pairedbin.rightExons  # Right exons of paired bin
        count = pairedbin.count  # Number of reads found

        for i in left:
            for j in right:
                if i == j:
                    right.remove(j)

        # If right exon list is now empty, no enu
        if len(right) == 0:
            new_bin = BinT(exons=left, count=count)
            if (len(new_bin.exons) > 2):
                all_pairedbins.append(new_bin)
            continue

        # If first right exon smaller than last left exon, no enu
        if right[0] < left[-1]:
            continue

        # The start node and end node for the enumeration has to be determined
        start_exon = left[len(left) - 1]
        end_exon = right[0]
        #print("start: ")
        #print(start_exon)
        for edgeKey, edgeValue in graph.edges.items():
            if edgeValue['type'] == 'Exon':
                if edgeValue['exon'] == start_exon:
                    start_node = edgeKey[0]
                elif edgeValue['exon'] == end_exon:
                    end_node = edgeKey[0]
            if start_node != "" and end_node != "":
                break

        # If start exon or end exon is not found, continue with next paired bin
        if start_node == "" or end_node == "":
            continue

        # Path enumeration between left and right exons is carried out to find all possible connections
        new_pairedbins = path_enumeration.enumeration_bins2(graph, new_pairedbins, start_node, [start_node], [],
                                                            multibins, end_node)

        # Bins have to be completed by adding the missing start and end exons
        for i in range(len(new_pairedbins)):
            if len(left) > 1:
                exons = left[0:-1] + new_pairedbins[i] + right
            else:
                exons = new_pairedbins[i] + right
            # Bin is added to the bin list
            new_bin = BinT(exons=exons, count=count)
            if (len(new_bin.exons) > 2):
                all_pairedbins.append(new_bin)
    return all_pairedbins





