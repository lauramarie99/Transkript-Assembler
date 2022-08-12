# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy

def writeGStar (graph:dict):
    graphStar = nx.DiGraph() # Define new Graph
    graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
    
    s0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set s0 label
    graphStar.add_node(s0Label) # Add s0
    graphStar.add_edge(s0Label, '0', weight = 0, capacity=float('inf')) # add s0 -> s
    
    t0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set t0 label
    graphStar.add_node(t0Label) # Add t0
    graphStar.add_edge(t0Label, '1', weight = 0, capacity=float('inf')) # add t -> t0

    graphStar.add_edge(t0Label, s0Label, weight = 0, capacity=float('inf')) # add t0 -> s0    
    sStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set sStarLabel
    graphStar.add_node(sStarLabel) # add s*
    
    tStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set tStarLabek
    graphStar.add_node(tStarLabel) # add t*

    graphStar.add_edges_from(graph.edges.keys(), capacity = float('inf')) # Add Edges from original graph with infinite capacity and cost function
    for edgeKey, edgeValue in graph.edges.items():
        graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = float(edgeValue['counts']['c'])) # Add backward edges with capacity (counts) of original forward edges

    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, tStarLabel, weight = 0,capacity = d)
        elif d<0:
            graphStar.add_edge(sStarLabel, node, capacity = -d, weight = 0)
    nx.draw(graphStar, with_labels=True)