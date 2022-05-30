#Work in progress. Script written so it can be excecuted directly. 
import Parsing
import networkx as nx
import FinalEnu
import new_path
unknown_bins = []
with open(
        "/home/tezie/info2/Test4.graph") as f:  # Test.graph consists of one single gene entry
    fileEndReached = False
    f.readline()
    while not fileEndReached:
        f.readline()
        Chromosome, Strand, Exons = Parsing.parse_meta(f)
        Bins = Parsing.parse_bins(f)
        multi_bins = FinalEnu.get_multibins(Bins)
        #multi_bins = Parsing.parse_bins(f)
        #print("Bins:", Bins)
        PairedBins = Parsing.parse_pairs(f)
        # BUILD GRAPHS
        G_full = nx.DiGraph()
        fileEndReached, skip = Parsing.parse_graph(f, G_full, Exons)  # Full Graph

        if not fileEndReached and not skip:
            G_clean = nx.DiGraph()
            fileEndReached, _ = Parsing.parse_graph(f, G_clean, Exons)  # Cleaned Graph
            # nx.draw_networkx(G_full, with_labels=True, arrowsize=12)
            # plt.show()
        print(len(PairedBins))
        print(PairedBins)
        #remove useless info i.e. if bins share same exons
        for bin in PairedBins:
            for i in bin.leftExons:
                for j in bin.rightExons:
                    if i == j:
                        bin.rightExons.remove(j)
        print(PairedBins)
        #create multibins where possible
        for bin in PairedBins:
            if len(bin.leftExons) != 0 and len(bin.rightExons) != 0:
                if ((bin.leftExons[-1]+1) == bin.rightExons[0]):
                    for i in bin.rightExons:
                        bin.leftExons.append(i)

                # get the unknown bins (i.e. all possibilities in the gaps)
                if((bin.leftExons[-1]+1) < bin.rightExons[0]):
                    transkripts_bins = []
                    n = bin.leftExons[-1]
                    d = bin.rightExons[0]
                    print(d)
                    transkripts_bins = FinalEnu.enumeration_bins2(G_clean, transkripts_bins,"0", ["0"], [],
                                                                           multi_bins, str(d))
                    #print(transkripts_bins)
                    unknown_bins.append(bin.leftExons)
                    unknown_bins.append((bin.rightExons))
                    #bin.leftExons.append(bin.rightExons)
        print(PairedBins)
        print(unknown_bins)
        print(transkripts_bins)

f.close()



