# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy
import math
import sys, os
import numpy as np
import scipy as sp

def writeGStar(graph:dict, costIndex:int, maxAdditionalEdgeCount:int, geneCounter):
    if costIndex<=2:
        return writeGStarLinear(graph, costIndex, geneCounter)
    else:
        return writeGStarQuadratic(graph, costIndex, maxAdditionalEdgeCount, geneCounter)

def writeGStarLinear(graph:dict, costIndex:int, geneCounter):
    graphStar = nx.DiGraph() # Define new Graph
    graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
    graphStar.add_node('s0') # Add s0
    graphStar.add_node('t0') # Add t0
    graphStar.add_node('s*') # add s*
    graphStar.add_node('t*') # add t*
    
    graphStar.add_edge('s0', '0', weight = 0, capacity=float('inf')) # add edge s0 -> s
    graphStar.add_edge('1', 't0', weight = 0, capacity=float('inf')) # add edge t -> t0
    graphStar.add_edge('t0', 's0', weight = 0, capacity=float('inf')) # add egde t0 -> s0    

    # Define scalingFactor to prevent floating point numbers in 1/cov(u,v) or 1/sqrt{cov(u,v)} for minCostFlow Calculation
    if costIndex == 0:
        scalingFactor = 1
    else:
        scalingFactor = 10e7
    # Add Coverage to Helper Edges
    for edgeKey, edgeValue in graph.edges.items():
        if edgeValue['type'] == 'Helper':
            if edgeKey[0] =='0':
                endNode = edgeKey[1]
                if (len(graph.out_edges(endNode, data=True)))>1: # Check for problems in the graph (i. e. start node of an exon has more than one successor)
                    print('More than one following edge')
                for edge in graph.out_edges(endNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c'] # Set Coverage/Count of source edges to coverage/count of succeeding edges
            elif edgeKey[1] =='1':
                startNode = edgeKey[0]
                if (len(graph.in_edges(startNode, data=True)))>1: # Check for problems in the graph (i. e. end node of an exon has more than one predecessor)
                    print('More than one following edge')
                for edge in graph.in_edges(startNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c'] # Set Coverage/Count of drain Edges to coverage/count of preceding edges

    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, 't*', capacity = d, weight = 0)
        elif d<0:
            graphStar.add_edge('s*', node, capacity = -d, weight = 0)

    # Now add inner edges to the graph
    for edgeKey, edgeValue in graph.edges.items():
        # Read edge Properties
        coverage = edgeValue['counts']['c']
        type = edgeValue['type']
        length = edgeValue['length']
        graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = float('inf'), weight=int(scalingFactor*(costFunction(1, 0, coverage, costIndex, length, type)))) # Add Forward edge with capacity infinite and weight = costFunction
        graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage, weight=int(scalingFactor*(costFunction(1, 0, coverage, costIndex, length, type)))) # Add backward edges with capacity (counts) of original forward edges and costFunction
    
    # Calculate max_flow at min cost
    flowDict = nx.max_flow_min_cost(graphStar, 's*', 't*', 'capacity', 'weight')
    # Correct cov(u,v) on original graph
    for edgeKey, edgeValue in graph.edges.items():
        edgeValue['counts']['c'] = edgeValue['counts']['c'] + flowDict[edgeKey[0]][edgeKey[1]] - flowDict[edgeKey[1]][edgeKey[0]]
    
    # Calculate total flow
    flow = 0 
    for edge in graph.out_edges('0', data=True):
        flow = flow + edge[2]['counts']['c']
    return(graphStar, graph, flow)    

def writeGStarQuadratic(graph:dict, costIndex:int, maxAdditionalEdgeCount, geneCounter):
    graphStar = nx.MultiDiGraph() # Define new MultiDiGraph
    graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
    graphStar.add_node('s0') # Add s0    
    graphStar.add_node('t0') # Add t0
    graphStar.add_node('s*') # add s*
    graphStar.add_node('t*') # add t*

    graphStar.add_edge('s0', '0', capacity=float('inf'), weight = 0) # add edge s0 -> s
    graphStar.add_edge('1', 't0', capacity=float('inf'), weight = 0) # add egde t -> t0
    graphStar.add_edge('t0', 's0', capacity=float('inf'),weight = 0) # add egde t0 -> s0    
    
    #Add Coverage to Helper Edges
    for edgeKey, edgeValue in graph.edges.items():
        if edgeValue['type'] == 'Helper':
            if edgeKey[0] =='0':
                endNode = edgeKey[1]
                if (len(graph.out_edges(endNode, data=True)))>1: # Check for problems in the graph (i. e. start node of an exon has more than one successor)
                    print('More than one following edge')
                for edge in graph.out_edges(endNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c']
            elif edgeKey[1] =='1':
                startNode = edgeKey[0]
                if (len(graph.in_edges(startNode, data=True)))>1: # Check for problems in the graph (i. e. end node of an exon has more than one predecessor)
                    print('More than one following edge')
                for edge in graph.in_edges(startNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c'] # Set Coverage/Count of drain Edges to coverage/count of preceding edges
    
    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, 't*', capacity = d, weight = 0)
        elif d<0:
            graphStar.add_edge('s*', node, capacity = -d, weight = 0)
    
    # Define Demand for s* and t*
    sourceDemand = sum(int(sourceEdge[2]['capacity']) for sourceEdge in graphStar.out_edges('s*', data=True))
    drainDemand = sum(int(drainEdge[2]['capacity']) for drainEdge in graphStar.in_edges('t*', data=True))
    
    
    # Add additional edges for every edge according to coverage
    
    # CostFunction 3: f(x) = x^2 modelled as ((i+x)*(i+x) - i*i)/x; 
    # CostFunction 4: f(x) = x^2/cov(u,v) modelled as ((i+x)*(i+x) - i*i)/(x*cov(u,v))
    # CostFunction 5: f(x( = x^2 * length(u,v)/cov(u,v) modelled as ((i+x)*(i+x) - i*i)*length(u,v)/(x*cov(u,v))

    edgeCounter = 0 
    if costIndex <6:
        for edgeKey, edgeValue in graph.edges.items():    
            # Read edge properties
            coverage = int(edgeValue['counts']['c']) # Assign Coverage of the particular edge to edge
            length = max(1, edgeValue['length']) 
            type = edgeValue['type'] 
            # Define scaling factor according to the used costfuntion
            if costIndex == 3: 
                scalingFactor = 1 
            else:
                scalingFactor = 10e7 
            maxAdditionalEdgeCount = min(100, int(math.ceil(25000/(len(graph.edges())))))
            # Calculate necessary stepSize (y and x) for each forward and backward edge to achieve at maximum maxAdditionalEdgeCount edges 
            y = int(max(1,math.floor(sourceDemand/maxAdditionalEdgeCount))) # Stepsize for ForwardEdges
            x = int(max(1, math.floor(coverage/maxAdditionalEdgeCount))) # Stepsize for Backwardedges
            for i in range(0,(min(sourceDemand, maxAdditionalEdgeCount))):
                edgeCounter +=1
                graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = y, weight=int(scalingFactor*(costFunction(i, y, coverage, costIndex, length, type)))) # Add Forward edge with capacity y and weight = costFunction
            for i in range(0, min(coverage, maxAdditionalEdgeCount)):
                edgeCounter +=1
                graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = x, weight=int(scalingFactor*(costFunction(i, x, coverage, costIndex, length, type)))) # Add backward edges with capacity x and costFunction    


            # Alternative soultion with correct modelling

            # y = int(max(1,math.ceil(sourceDemand/maxAdditionalEdgeCount))) # Stepsize for ForwardEdges
            # x = int(max(1, math.ceil(coverage/maxAdditionalEdgeCount))) # Stepsize for Backwardedges
            # for i in range(0,(min(sourceDemand, maxAdditionalEdgeCount))):
            #     if i*y +y <= sourceDemand:
            #         edgeCounter +=1
            #         graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = y, weight=int(scalingFactor*(costFunction(i, y, coverage, costIndex, length, type)))) # Add Forward edge with capacity y and weight = costFunction
            #         if i*y +y == sourceDemand:
            #             break
            #     else:
            #         graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = sourceDemand - y*i, weight=int(scalingFactor*(costFunction(i, sourceDemand - y*i, coverage, costIndex, length, type)))) # Add Forward edge with capacity y and weight = costFunction
            #         break
            # for i in range(0, min(coverage, maxAdditionalEdgeCount)):
            #     if i*x + x <= coverage:
            #         edgeCounter +=1
            #         graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = x, weight=int(scalingFactor*(costFunction(i, x, coverage, costIndex, length, type)))) # Add backward edges with capacity x and costFunction    
            #         if i*x + x == coverage:
            #             break
            #     else:
            #         graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage -x*i, weight=int(scalingFactor*(costFunction(i, coverage-x*i, coverage, costIndex, length, type)))) # Add backward edges with capacity x and costFunction    
            #         break

    # Not recommended Costfunctions because of heavy computation time
        # Costfunction 6: f(x) = x^2 (modeled as ((i+1)^2-i^2) 
        # Costfunction 7: f(x) = x^2 (modelled as x(x+1)/2) 
        # Costfunction 8: f(x) = x^2/cov(u,v) modelled as x(x+1)/2)/cov(u,v) 

    else: 
        scalingFactor=1
        if costIndex == 8:
            scalingFactor=10e7  # Cost-Function: 3: i; 4:((i+1)*(i+1) - i*i) -> requires very long computation time, since no stepSize is included
        for edgeKey, edgeValue in graph.edges.items():    
            coverage = int(edgeValue['counts']['c']) # Assign Coverage of the particular edge to edge
            type = edgeValue['type'] # Assign edge type to the variable type
            length = max(1, edgeValue['length']) # Assign length of the edge to the variable "length", Helper edges will receive the length 0
            for i in range(sourceDemand):
                graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = 1, weight=int(scalingFactor*(costFunction(i, 1, coverage, costIndex, length, type)))) # Add Forward edge with capacity 1 and weight = costFunction
            for i in range(coverage):
                graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = 1, weight=int(scalingFactor*(costFunction(i, 1, coverage, costIndex, length, type)))) # Add backward edges with capacity 1 and costFunction    

    # Write DemandDictionary
    nodeDemand = {}
    for nodeKey in graphStar.nodes.keys():
        if nodeKey == 't*':
            nodeDemand [nodeKey] = drainDemand
        elif nodeKey == 's*':
            nodeDemand [nodeKey] = -sourceDemand
        else:
            nodeDemand[nodeKey] = 0
    
    # Set Demand for each node in Off-Set Network
    
    nx.set_node_attributes(graphStar, nodeDemand, name='demand')
    
    # Calculate min_cost_flow
    flowDict = nx.min_cost_flow(graphStar, 'demand', 'capacity', 'weight')
    # Correct on cov(u,v) on original graph
    for edgeKey, edgeValue in graph.edges.items():
        edgeFlowForward = sum (flowDict[edgeKey[0]][edgeKey[1]][i] for i in range(len(flowDict[edgeKey[0]][edgeKey[1]]))) # Sum flow on all forward edges
        edgeFlowBackward = sum (flowDict[edgeKey[1]][edgeKey[0]][i] for i in range(len(flowDict[edgeKey[1]] [edgeKey[0]]))) # Sum flow on all backward edges
        edgeValue['counts']['c'] = edgeValue['counts']['c'] + edgeFlowForward - edgeFlowBackward # calculate flow in the original graph
    
    # Calculate total flow 
    flow = 0 
    #totalFlow = sum(edge[2]['counts']['c'] for edge in graph.out_edges('0', data=True))
    for edge in graph.out_edges('0', data=True):
        flow = flow + edge[2]['counts']['c']
    return(graphStar, graph, flow)    

# Specify costFunctions
def costFunction (i, x: int, coverage:int, costIndex, length:int, type:str):
    if costIndex == 0:
        return 1
    elif costIndex == 1:
        return 1/coverage
    elif costIndex == 2: 
        return 1/math.sqrt(coverage)
    elif costIndex == 3:
        return ((i+x)*(i+x) - i*i)/x
    elif costIndex == 4:
        return ((i+x)*(i+x) - i*i)/(coverage*x)
    elif costIndex == 5:
        return ((i+x)*(i+x) - i*i)*length/(coverage*x)
    elif costIndex == 6:
        return ((i+1)*(i+1) - i*i)
    elif costIndex == 7:
        return i
    elif costIndex == 8:
        return i/coverage

def flowDecompositionWithTranscriptlist(graph:dict, transcriptsCopy:list, decomposition_option:str, flow): 
    # define new list for optimized transcripts
    optimizedTranscripts = []
    # Compute FlowDecomposition by removing the flow of the longest path from the transcript list
    if decomposition_option == 'longestPath': 
        while (flow!=0 and len(transcriptsCopy)>0): 
            longestPath = transcriptsCopy[transcriptsCopy.index(max(transcriptsCopy, key=len))] # get longest Path from transcriptsCopy List
            minFlow = min(graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] for i in range(len(longestPath)-1)) # Define minimal flow on this path = pathFlow
            if flow - minFlow >= 0 and minFlow>0: # If removal of this flow results in non-negative flow
                for i in range(len(longestPath)-1): 
                    graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] = graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] - minFlow # reduce flow of every edge by this flow
                flow = flow - minFlow # reduce TotalFlow by pathFlow
                optimizedTranscripts.append((longestPath, minFlow)) # Add path to optimized TranscriptList
            transcriptsCopy.remove(longestPath) # Remove path from transcriptsCopy
    
    # Compute FlowDecomposition by removing the flow of the path with maximumFlow from the transcript list
    elif decomposition_option == 'maximumFlow':
        while (flow!=0 and len(transcriptsCopy)>0):
            maxFlow = 0 # set maxFlow to 0
            maxFlowPath = [] # define new List with maxFlowPath
            for transcript in transcriptsCopy:
                # Read the minimum flow for every path of trancsriptsCopy 
                minFlow = min(graph.edges[transcript[i], transcript[i+1]]['counts']['c'] for i in range(len(transcript)-1)) 
                if minFlow > maxFlow: # if it's bigger than the already found flow -> overwrite it
                    maxFlow = minFlow 
                    maxFlowPath = transcript
            if ((flow - maxFlow) >= 0) and flow>0: # if the remaining total flow is non-negative remove this flow and add the path to optimized transcripts  
                for i in range(len(maxFlowPath)-1):
                    graph.edges[maxFlowPath[i], maxFlowPath[i+1]]['counts']['c'] = graph.edges[maxFlowPath[i], maxFlowPath[i+1]]['counts']['c'] - maxFlow
                flow = flow - maxFlow
                optimizedTranscripts.append((maxFlowPath, maxFlow))
            transcriptsCopy.remove(maxFlowPath) # remove path from the transcriptList in any case 
    else:
        print('Error, please specify decomposition option.')
        os.exit()
    return optimizedTranscripts, flow

def flowDecompositionDP (graph: dict, decomposition_option:str, flow:int, geneCounter):
    optimizedTranscripts = []
    if decomposition_option == 'longestPath':
        while(flow!=0):
            # Forward
            longestPathDict = {}
            queue = [] # Initialize Queue
            queue.append('0') # Append first Item
            longestPathDict['0'] = (None, 0) # Set length of source to successors to 1
            queue = [] # Initialize PriorityQueue
            for v in list(nx.topological_sort(graph)):
                length=-1
                for u in list(graph.predecessors(v)):
                    if u in longestPathDict.keys() and length <longestPathDict[u][1]:
                        length = longestPathDict[u][1]
                        longestPathDict[v] = (u, longestPathDict[u][1]+1)
            # Backtracking
            transcript = []
            length = longestPathDict['1'][1] 
            index = longestPathDict['1'][0]
            transcript.append('1')
            while(length>0):
                transcript.append(index)
                length = longestPathDict[index][1]
                index = longestPathDict[index][0]
            transcript.reverse()

            # Eliminate edges with minimal flow
            minFlow = min(graph.edges[transcript[i], transcript[i+1]]['counts']['c'] for i in range(len(transcript)-1)) # Define minFlow
            for i in range(len(transcript)-1):
                graph.edges[transcript[i], transcript[i+1]]['counts']['c'] = graph.edges[transcript[i], transcript[i+1]]['counts']['c'] - minFlow 
                if graph.edges[transcript[i], transcript[i+1]]['counts']['c'] == 0:
                    graph.remove_edge(transcript[i], transcript[i+1])
            flow = flow - minFlow
            optimizedTranscripts.append((transcript, minFlow))

    elif decomposition_option == 'maximumFlow':
        topologicalOrder = list(nx.topological_sort(graph))
        while(flow!=0):        
            #Forward 
            maxFlowDict = {}
            for key in graph.edges.keys():
                maxFlowDict[key] = graph.edges[key]['counts']['c']
            for v in topologicalOrder:
                if v=='0':
                    continue
                maxPreFlow = max(maxFlowDict[(x,v)] for x in graph.predecessors(v))
                for u in list(graph.adj[v]):
                    maxFlowDict[(v,u)] = min(maxFlowDict[(v,u)], maxPreFlow)
                if v=='1':
                    break
            if max(maxFlowDict[(x,'1')] for x in graph.predecessors('1')) == 0:
                break

            # Backtracking
            else:
                transcript = []
                transcript.append('1')
                minMaxFlow = 0
                for x in graph.predecessors('1'):
                    if maxFlowDict[(x,'1')] > minMaxFlow:
                        maxIndex = x
                        minMaxFlow = maxFlowDict[(x,'1')]
                transcript.append(str(maxIndex))
                while(maxIndex!='0'):
                    u = maxIndex
                    maxFlow = 0
                    for x in graph.predecessors(u):
                        if maxFlowDict[(x,u)] > maxFlow:
                            maxFlow = maxFlowDict[(x, u)]
                            maxIndex = x
                    transcript.append(str(maxIndex))
                transcript.reverse()
                optimizedTranscripts.append((transcript, minMaxFlow))
                
                # Eliminate edges with minimal flow
                for i in range(len(transcript)-1):
                    graph.edges[transcript[i], transcript[i+1]]['counts']['c'] = graph.edges[transcript[i], transcript[i+1]]['counts']['c'] - minMaxFlow
                flow = flow - minMaxFlow  
    return (optimizedTranscripts, flow)

