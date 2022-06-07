#get this shit cleaned up by tomorrow....


# IMPORT
from collections import namedtuple
from pydoc import pathdirs
import Parsing as parse_graph_new
import new_path as path_enumeration
import MyEnu
import networkx as nx
from matplotlib import pyplot as plt
import time
from matplotlib import pyplot as plt
start_time = time.time()
# MAIN

total_trans = 0
with open(
        "/home/tezie/info2/real.graph") as f:  # Test.graph consists of one single gene entry
    fileEndReached = False
    f.readline()
    while not fileEndReached:
        f.readline()
        Chromosome, Strand, Exons = parse_graph_new.parse_meta(f)
        Bins = parse_graph_new.parse_bins(f)
        #print("Bins:", Bins)
        PairedBins = parse_graph_new.parse_pairs(f)

        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = parse_graph_new.parse_graph(f, G_full, Exons)  # Full Graph

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()
            fileEndReached, _ = parse_graph_new.parse_graph(f, G_clean, Exons)  # Cleaned Graph
            #nx.draw_networkx(G_full, with_labels=True, arrowsize=12)
            #plt.show()

        # PATH ENUMERATION OF GRAPH

        transkripts_bins1 = []  # All transkripts obtained by enumeration function 1
        transkripts_bins2 = []  # All transkripts obtained by enumeration function 2
        transkripts_bins3 = []
        transkripts_bins4 = []
        multi_bins = path_enumeration.get_multibins(Bins)  # Filter all bins with more than two exons

        #transkripts_bins1 = path_enumeration.enumeration_bins(G_clean, transkripts_bins1, "0", ["0"], multi_bins,
                                                           #   multi_bins, "1")  # Function 1
        #print("Transkripts for function 1:", len(transkripts_bins1), transkripts_bins1)

        transkripts_bins2 = path_enumeration.enumeration_bins2(G_full, transkripts_bins2, "0", ["0"], [], multi_bins,
                                                               "1")
        total_trans = total_trans + len(transkripts_bins2) # Function 2
        # print("Transkripts for function 2:", len(transkripts_bins2), transkripts_bins2)
        #
        # transkripts_bins3 = MyEnu.EnuBinConstraint(G_clean, transkripts_bins3, ["0"], "0", [], multi_bins,"1"
        #                                                        )
        # print("Transkripts for function 3:", len(transkripts_bins3), transkripts_bins3)
        # # The results should be the same!
        # transkripts_bins4 = path_enumeration.enumeration(G_clean, transkripts_bins4, "0", ["0"],"1")
        # print("Transkripts for function 4:", len(transkripts_bins4), transkripts_bins4)

f.close()

print(total_trans)
print("--- %s seconds ---" % (time.time() - start_time))
