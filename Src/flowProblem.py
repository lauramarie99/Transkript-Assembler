# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy
import math
import sys, os
import numpy as np
import scipy as sp

def writeGStar (graph:dict, costIndex):
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
        # Add Coverage to Helper Edges
        if edgeValue['type'] == 'Helper':
            if edgeKey[0] =='0':
                endNode = edgeKey[1]
                for edge in graph.out_edges(endNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c']
            elif edgeKey[1] =='1':
                startNode = edgeKey[0]
                for edge in graph.in_edges(startNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c']
        coverage = edgeValue['counts']['c']
        graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = float('inf'), weight=costFunction(1, float('inf'), costIndex, 0)) # Add Forward edge with capacity infinite and weight = costFunction
        if costIndex <= 2:
            graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage, weight=costFunction(1, coverage, costIndex, 0)) # Add backward edges with capacity (counts) of original forward edges and costFunction
        else: 
            for i in range(1,coverage+1):
                graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = 1, weight=costFunction(i, 1, costIndex, 0)) # Add backward edges with capacity (counts) of original forward edges and costFunction
    print(len(graphStar.edges()))

    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, 't*', weight = 0,capacity = d)
        elif d<0:
            graphStar.add_edge('s*', node, weight = 0, capacity = -d)

    # Calculate maximum Flow at minimal Costs 
    flowDict = nx.max_flow_min_cost(graphStar, 's*', 't*', 'capacity', 'weight')

    # Correct on cov(u,v) on original graph
    for edgeKey, edgeValue in graph.edges.items():
        edgeValue['counts']['c'] = edgeValue['counts']['c'] + flowDict[edgeKey[0]][edgeKey[1]] - flowDict[edgeKey[1]][edgeKey[0]]
    flow = 0 
    for edge in graph.out_edges('0', data=True):
        flow = flow + edge[2]['counts']['c']
    return(graphStar, graph, flow)    

def costFunction (i, coverage:int, costIndex, length):     
    if costIndex == 0:
        return 1
    elif costIndex == 1:
        return 1/coverage
    elif costIndex == 2: 
        return 1/math.sqrt(coverage)
    elif costIndex == 3:
        return i/coverage
    elif costIndex == 4:
        return (abs(length)-1)/coverage 
    else:
        print('Please enter valid CostFunctionIndex, for specifications see main.py --help. The programm will now be terminated.')
        sys.exit()

def flowDecompositionWithTranscriptlist(graph:dict, transcripts:list, decomposition_option:str, flow):
    transcriptsCopy = deepcopy(transcripts)
    #transcriptsCopy = list(nx.all_simple_paths(graph, '0', '1'))
    optimizedTranscripts = []
    if decomposition_option == 'longestPath':
        while (flow>0 and len(transcriptsCopy)>0):
            longestPath = transcriptsCopy[transcriptsCopy.index(max(transcriptsCopy, key=len))]
            minFlow = min(graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] for i in range(len(longestPath)-1))
            if flow - minFlow >= 0:
                for i in range(len(longestPath)-1):
                    graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] = graph.edges[longestPath[i], longestPath[i+1]]['counts']['c'] - minFlow 
                flow = flow - minFlow
                optimizedTranscripts.append(longestPath) # Add path to optimized TranscriptList
            transcriptsCopy.remove(longestPath) # Remove path from transcriptsCopy
        print('flowDecompositionWithTranscriptlist and longestPath ' + str(flow))
    elif decomposition_option == 'maximumFlow':
        while (flow>0 and len(transcriptsCopy)>0):
            maxFlow = 0
            maxFlowPath = []
            flowCounter = 0
            for transcript in transcriptsCopy:
                minFlow = min(graph.edges[transcript[i], transcript[i+1]]['counts']['c'] for i in range(len(transcript)-1))
                if minFlow > maxFlow:
                    maxFlow = minFlow
                    maxFlowPath = transcript
                    flowCounter = 1
            if flowCounter == 0:
                break
            if ((flow - maxFlow) >= 0) :
                for i in range(len(maxFlowPath)-1):
                    graph.edges[maxFlowPath[i], maxFlowPath[i+1]]['counts']['c'] = graph.edges[maxFlowPath[i], maxFlowPath[i+1]]['counts']['c'] - maxFlow 
                flow = flow - maxFlow
                optimizedTranscripts.append(maxFlowPath)
            transcriptsCopy.remove(maxFlowPath)
        print('flowDecompositionWithTranscriptlist and maxFlow ' + str(flow))
    else:
        print('Error, please specify decomposition option.')
        os.exit()
    return optimizedTranscripts

def flowDecompositionDP (graph: dict, decomposition_option:str, flow):
    optimizedTranscripts = []
    if decomposition_option == 'longestPath':
        while(flow>0):
            # Forward
            n = max(int(node) for node in list(graph.nodes()))+1 # Define length of matrix
            DP = np.zeros((n,n), dtype=int) # Initialize Matrix
            queue = [] # Initialize Priority Queue
            queue.append('0') # Append first Item
            v = queue.pop(0) # Pop item an read first item
            succ = list(graph.adj[v]) 
            for u in succ: 
                DP[int(v)][int(u)] = 1 
                queue.append(u)
            while(len(queue)>0):
                v = queue.pop(0)
                succ = list(graph.adj[v])
                for u in succ:
                    DP[int(v)][int(u)] = max(DP[int(x)][int(v)] for x in graph.predecessors(v)) + 1 
                    queue.append(u)

            # Backtracking
            transcript = []
            maxIndex = np.argmax(DP, axis=0)[1]
            length = DP[maxIndex][1]
            transcript.append('1')
            while(length>0):
                v = maxIndex
                transcript.append(str(maxIndex))
                maxIndex = np.argmax(DP, axis=0)[v]
                length = DP[maxIndex][v]
            transcript.reverse()
            optimizedTranscripts.append(transcript)
            
            # Eliminate edges with minimal flow
            minFlow = min(graph.edges[transcript[i], transcript[i+1]]['counts']['c'] for i in range(len(transcript)-1)) # Define minFlow
            for i in range(len(transcript)-1):
                graph.edges[transcript[i], transcript[i+1]]['counts']['c'] = graph.edges[transcript[i], transcript[i+1]]['counts']['c'] - minFlow 
                if graph.edges[transcript[i], transcript[i+1]]['counts']['c'] == 0:
                    graph.remove_edge(transcript[i], transcript[i+1])
            flow = flow - minFlow
        print('flowDecompositionWithDP and longestPath ' + str(flow))

    elif decomposition_option == 'maximumFlow':
        while(flow>0):        
            # Forward 
            n = max(int(node) for node in list(graph.nodes()))+1
            DP = np.zeros((n,n), dtype=int) # Initialize DP - Matrix
            queue = [] # Initialize PriorityQueue
            for u in list(graph.adj['0']):
                DP[0][int(u)]= graph.edges['0', u]['counts']['c']
                queue.append(u)
            while(len(queue)>0):
                v = queue.pop(0)
                maxPreFlow = max(DP[int(x)][int(v)] for x in graph.predecessors(v))
                for u in list(graph.adj[v]):
                    DP[int(v)][int(u)] = min(graph.edges[v,u]['counts']['c'], maxPreFlow)
                    queue.append(u)

            # Backtracking
            transcript = []
            transcript.append('1')
            maxFlow = max(DP[int(x)][1] for x in graph.predecessors('1'))
            maxIndex = np.argmax(DP, axis=0)[1]
            transcript.append(str(maxIndex))
            while(maxIndex!=0):
                maxIndex = np.argmax(DP, axis=0)[maxIndex]
                transcript.append(str(maxIndex))
            transcript.reverse()
            optimizedTranscripts.append(transcript)
            
            # Eliminate edges with minimal flow
            for i in range(len(transcript)-1):
                graph.edges[transcript[i], transcript[i+1]]['counts']['c'] = graph.edges[transcript[i], transcript[i+1]]['counts']['c'] - maxFlow 
                if graph.edges[transcript[i], transcript[i+1]]['counts']['c'] == 0:
                    graph.remove_edge(transcript[i], transcript[i+1])
            flow = flow - maxFlow
        print('flowDecompositionWithDP and residualFlow = ' + str(flow))

    # Second Option with One matrix 
        # n = max(int(node) for node in list(graph.nodes()))+1
        # DP = np.zeros((n,n), dtype=int) # Initialize DP - Matrix
        # queue = [] # Initialize PriorityQueue
        # queue.append('0') # Add source to Priority Queue 
        # while(len(queue)>0): # Fill DP-Matrix
        #     v = queue.pop(0) 
        #     for u in list(graph.adj[v]):
        #         queue.append(u)
        #         DP[int(v)][int(u)] = graph.edges[v,u]['counts']['c']
        
        # # Backtracing
        # while(flow>0):
        #     transcript = []
        #     transcript.append('1')
        #     maxFlow = max(DP[int(x)][1] for x in graph.predecessors('1'))
        #     maxIndex = np.argmax(DP, axis=0)[1]
        #     transcript.append(str(maxIndex))
        #     while(maxIndex!=0):
        #         edgeFlow =max(DP[int(x)][maxIndex] for x in graph.predecessors(str(maxIndex)))
        #         if edgeFlow < maxFlow:
        #             maxFlow = edgeFlow
        #         maxIndex = np.argmax(DP, axis=0)[maxIndex]
        #         transcript.append(str(maxIndex))
        #     transcript.reverse()
        #     optimizedTranscripts.append(transcript)
            
        #     # Update DP-Matrix with reduced flow
        #     for i in range(len(transcript)-1):
        #         DP[int(transcript[i])][int(transcript[i+1])] = DP[int(transcript[i])][int(transcript[i+1])] - maxFlow
        #     flow = flow - maxFlow
    return optimizedTranscripts

