# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy
import math

def writeGStar (graph:dict):
    graphStar = nx.DiGraph() # Define new Graph
    graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
    
    #s0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set s0 label
    graphStar.add_node('s0') # Add s0
    graphStar.add_edge('s0', '0', weight = 0, capacity=float('inf')) # add s0 -> s
    
    #t0Label = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set t0 label
    graphStar.add_node('t0') # Add t0
    graphStar.add_edge('t0', '1', weight = 0, capacity=float('inf')) # add t -> t0

    graphStar.add_edge('t0', 's0', weight = 0, capacity=float('inf')) # add t0 -> s0    
    #sStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # set sStarLabel
    graphStar.add_node('s*') # add s*
    
    #tStarLabel = str(max(int(nodeNumber) for nodeNumber in graphStar.nodes.keys())+1) # Set tStarLabek
    graphStar.add_node('t*') # add t*

    # Add Edges from original graph with infinite capacity and cost function

    for edgeKey, edgeValue in graph.edges.items():
        coverage = edgeValue['counts']['c']
        graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = float('inf'), weight=costFunction(coverage)) # Add Forward edge with capacity infinite and weight = costFunction
        graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage, weight=costFunction(coverage)) # Add backward edges with capacity (counts) of original forward edges and costFunction
        
    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            x = 0 # Fluss -> zu definieren
            graphStar.add_edge(node, 't*', weight =0,capacity = d)
        elif d<0:
            x =0 # Fluss -> zu definieren
            graphStar.add_edge('s*', node, weight = 0, capacity = -d)
    #nx.draw(graphStar, with_labels=True)
    
    # Define Coverage
    coverage1 = []
    for edgeValue in graphStar.edges.values():
        coverage1.append(float(edgeValue['capacity']))

    numberEdges = len(graphStar.edges)
    edges = list(graphStar.edges.keys())

    # Introduce Model for Minimizing the flow 
    flowModel = gp.Model()

    # Add decision Variables
    flow = flowModel.addVars(numberEdges, vtype = GRB.CONTINUOUS)
    
    # Add SlackVariables
    y = flowModel.addVars(numberEdges, lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
    z = flowModel.addVars(numberEdges, vtype = GRB.CONTINUOUS)
    
    for j in range(len(flow)):
        flowModel.addConstr(y[j]==flow[j] - coverage1[j])
        flowModel.addConstr(z[j]>= y[j])
        flowModel.addConstr(z[j]>= -y[j])
    
    flowModel.setObjective(gp.quicksum(z[j] for j in range(len(flow))))
    flowModel.optimize()

    # Add Constraints
    for var in flowModel.getVars():
        print(var.varName)
        print(var.X)
    print(flow)
#     flowModel.optimize()
    return(graphStar)

def costFunction (coverage:int):     
    cost = 0
    return cost