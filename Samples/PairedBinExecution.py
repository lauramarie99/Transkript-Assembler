#uses the second PairedBinEnu script
import Parsing
import Paired_bin_Enu
import FinalEnu
import networkx as nx
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
        Chromosome, Strand, Exons = Parsing.parse_meta(f)
        Bins = Parsing.parse_bins(f)

        PairedBins = Parsing.parse_pairs(f)

        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = Parsing.parse_graph(f, G_full, Exons)  # Full Graph

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()
            fileEndReached, _ = Parsing.parse_graph(f, G_clean, Exons)  # Cleaned Graph
            #nx.draw_networkx(G_clean, with_labels=True, arrowsize=12)
            #plt.show()

        # PATH ENUMERATION OF GRAPH

        multibins = FinalEnu.get_multibins(Bins)
        pairedbins = Paired_bin_Enu.get_pairedbins(G_clean, PairedBins, multibins)
        transkripts = FinalEnu.enumeration_bins2(G_clean, [], "0", ["0"], [], (pairedbins + multibins), "1")
        total_trans = len(transkripts) + total_trans
f.close()
print(total_trans)
print("--- %s seconds ---" % (time.time() - start_time))
