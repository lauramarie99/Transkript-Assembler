# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy
import math

def writeGStar (graph:dict, cost_index):
    graphStar = nx.DiGraph() # Define new Graph
    graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
    
    #s0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set s0 label
    graphStar.add_node('s0') # Add s0
    graphStar.add_edge('s0', '0', weight = 0, capacity=float('inf')) # add s0 -> s
    
    #t0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set t0 label
    graphStar.add_node('t0') # Add t0
    graphStar.add_edge('1', 't0', weight = 0, capacity=float('inf')) # add t -> t0
    graphStar.add_edge('t0', 's0', weight = 0, capacity=float('inf')) # add t0 -> s0    
    #sStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set sStarLabel
    graphStar.add_node('s*') # add s*
    
    #tStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set tStarLabek
    graphStar.add_node('t*') # add t*

    # Add Edges from original graph with infinite capacity and cost function

    for edgeKey, edgeValue in graph.edges.items():
        #Add Coverage to Helper Edges
        if edgeValue['type'] == 'Helper':
            if edgeKey[0] =='0':
                endNode = edgeKey[1]
                if (len(graph.out_edges(endNode, data=True)))>1:
                    print('More than one following edge')
                for edge in graph.out_edges(endNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c']
            elif edgeKey[1] =='1':
                startNode = edgeKey[0]
                if (len(graph.in_edges(startNode, data=True)))>1:
                    print('More than one following edge')
                for edge in graph.in_edges(startNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c']
        coverage = edgeValue['counts']['c']
        graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = float('inf'), weight=costFunction(float('inf'), cost_index, 0)) # Add Forward edge with capacity infinite and weight = costFunction
        graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage, weight=costFunction(coverage, cost_index, 0)) # Add backward edges with capacity (counts) of original forward edges and costFunction
        

    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, 't*', weight = 0,capacity = d)
        elif d<0:
            graphStar.add_edge('s*', node, weight = 0, capacity = -d)
    
    flowDict = nx.max_flow_min_cost(graphStar, 's*', 't*', 'capacity', 'weight')

    # Define Coverage
    coverage1 = {}
    for edgeKey,edgeValue in graphStar.edges.items():
        coverage1[edgeKey]= float(edgeValue['capacity'])

    # Correct on cov(u,v) on original graph
    for edgeKey, edgeValue in graph.edges.items():
        edgeValue['counts']['c'] = edgeValue['counts']['c'] + flowDict[edgeKey[0]][edgeKey[1]] - flowDict[edgeKey[1]][edgeKey[0]]
    
    # Flow Decomposition

    # Longest Path
    # for path in transcripts:
    #     print(transcripts)

    return(graphStar)    

def costFunction (coverage:int, cost_index, length):     
    if cost_index == 0:
        return 1
    elif cost_index == 1:
        print(coverage)
        print(1/coverage)
        return 1/coverage
    elif cost_index == 2: 
        return 1/math.sqrt(coverage)    