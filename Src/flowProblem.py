# Import packages
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from copy import deepcopy
import math
import sys, os
import numpy as np
import scipy as sp

def writeGStar(graph:dict, costIndex:int):
    if costIndex<=2:
        return writeGStarLinear(graph, costIndex)
    else:
        return writeGStarQuadratic(graph, costIndex)

def writeGStarLinear(graph:dict, costIndex:int):
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
    maxCoverage = max(edgeValue['counts']['c'] for edgeKey, edgeValue in graph.edges.items())

    # Define scalingFactor to prevent floating point numbers in 1/cov(u,v) or 1/sqrt{cov(u,v)}
    if costIndex == 0:
        scalingFactor = 1
    else:
        scalingFactor = maxCoverage * maxCoverage
    
    # Now add edges to the graph
    for edgeKey, edgeValue in graph.edges.items():
        #Add Coverage to Helper Edges
        if edgeValue['type'] == 'Helper':
            if edgeKey[0] =='0':
                endNode = edgeKey[1]
                if (len(graph.out_edges(endNode, data=True)))>1: # Check for problems in the graph (start node of an exon has more than one successor)
                    print('More than one following edge')
                for edge in graph.out_edges(endNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c'] # Set Coverage/Count of source edges to coverage/count of succeeding edges
            elif edgeKey[1] =='1':
                startNode = edgeKey[0]
                if (len(graph.in_edges(startNode, data=True)))>1: # Check for problems in the graph (end node of an exon has more than one predecessor)
                    print('More than one following edge')
                for edge in graph.in_edges(startNode, data=True):
                    edgeValue['counts']['c'] = edge[2]['counts']['c'] # Set Coverage/Count of drain Edges to coverage/count of preceding edges
        coverage = edgeValue['counts']['c'] # Define coverage variable and assign count to it
        type = edgeValue['type']
        graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = float('inf'), weight=int(scalingFactor*(costFunction(1, coverage, costIndex, 0, type)))) # Add Forward edge with capacity infinite and weight = costFunction
        graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = coverage, weight=int(scalingFactor*(costFunction(1, coverage, costIndex, 0, type)))) # Add backward edges with capacity (counts) of original forward edges and costFunction
        
    # Add additional edges from s* and to t*
    for node in graph.nodes.keys():
        d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
        if d>0:
            graphStar.add_edge(node, 't*', capacity = d, weight = 0)
        elif d<0:
            graphStar.add_edge('s*', node, capacity = -d, weight = 0)
    
    # Calculate max_flow at min cost
    print('Trying to calculate maximum flow at minimal costs')
    flowDict = nx.max_flow_min_cost(graphStar, 's*', 't*', 'capacity', 'weight')
    print('Finished calculating maximum flow at minimal costs')

    # Correct on cov(u,v) on original graph
    for edgeKey, edgeValue in graph.edges.items():
        edgeValue['counts']['c'] = edgeValue['counts']['c'] + flowDict[edgeKey[0]][edgeKey[1]] - flowDict[edgeKey[1]][edgeKey[0]]
    flow = 0 
    for edge in graph.out_edges('0', data=True):
        flow = flow + edge[2]['counts']['c']
    return(graphStar, graph, flow)    

def writeGStarQuadratic(graph:dict, costIndex:int):
    if sum(edgeValue['counts']['c'] for edgeKey, edgeValue in graph.edges.items())>50000:
        print('Number of additional edges is too high to calculate maximum Flow at minimum costs')
        return 0, 1, 2
    else:
        graphStar = nx.MultiDiGraph() # Define new Graph
        graphStar.add_nodes_from(graph.nodes.keys()) # Add all nodes from the original Graph
        graphStar.add_node('s0') # Add s0
        graphStar.add_edge('s0', '0', capacity=float('inf'), weight = 0) # add s0 -> s
        graphStar.add_node('t0') # Add t0
        graphStar.add_edge('1', 't0', capacity=float('inf'), weight = 0) # add t -> t0
        graphStar.add_edge('t0', 's0', capacity=float('inf'),weight = 0) # add t0 -> s0    
        graphStar.add_node('s*') # add s*
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
            coverage = int(edgeValue['counts']['c'])
            length = edgeValue['length']
            type = edgeValue['type']
            for i in range(1, coverage+1):
                # Add forward edges
                graphStar.add_edge(edgeKey[0], edgeKey[1], capacity = 1, weight=int(costFunction(i, coverage, costIndex, length, type))) # Add Forward edge with capacity infinite and weight = costFunction
                # Add Backward Edges
                graphStar.add_edge(edgeKey[1], edgeKey[0], capacity = 1, weight=int(costFunction(i, coverage, costIndex, length, type))) # Add backward edges with capacity (counts) of original forward edges and costFunction

        # Add additional edges from s* and to t*
        for node, nodeValue in graph.nodes.items():
            d = sum(int(edge[2]['counts']['c']) for edge in graph.out_edges(node, data=True)) - sum(int(edge[2]['counts']['c']) for edge in graph.in_edges(node, data=True))
            if d>0:
                graphStar.add_edge(node, 't*', capacity = d, weight = 0)
            elif d<0:
                graphStar.add_edge('s*', node, capacity = -d, weight = 0)
        
        # Define Demand for s* and t*
        sourceDemand = sum(int(sourceEdge[2]['capacity']) for sourceEdge in graphStar.out_edges('s*', data=True))
        drainDemand = sum(int(drainEdge[2]['capacity']) for drainEdge in graphStar.in_edges('t*', data=True))
        # Write DemandDictionary
        nodeDemand = {}
        for nodeKey, nodeValue in graphStar.nodes.items():
            if nodeKey == 't*':
                nodeDemand [nodeKey] = drainDemand
            elif nodeKey == 's*':
                nodeDemand [nodeKey] = -sourceDemand
            else:
                nodeDemand[nodeKey] = 0
        
        # Set Demand for each node in Off-Set Network
        nx.set_node_attributes(graphStar, nodeDemand, name='demand')
        
        # Calculate min_cost_flow
        print('Trying to calculate maximum Flow at minimal costs')
        flowDict = nx.min_cost_flow(graphStar, 'demand', 'capacity', 'weight')
        print('Finished calculating maximum Flow at minimal costs')
        
        # Correct on cov(u,v) on original graph
        for edgeKey, edgeValue in graph.edges.items():
            edgeFlowForward = sum (flowDict[edgeKey[0]][edgeKey[1]][i] for i in range(len(flowDict[edgeKey[0]][edgeKey[1]])))
            edgeFlowBackward = sum (flowDict[edgeKey[1]][edgeKey[0]][i] for i in range(len(flowDict[edgeKey[1]][edgeKey[0]])))  
            edgeValue['counts']['c'] = edgeValue['counts']['c'] + edgeFlowForward - edgeFlowBackward     
        flow = 0 
        flowDrain = 0
        for edge in graph.out_edges('0', data=True):
            flow = flow + edge[2]['counts']['c']
        for edge in graph.in_edges('1', data=True):
            flowDrain = flowDrain + edge[2]['counts']['c']
        return(graphStar, graph, flow)    
  

def costFunction (i, coverage:int, costIndex, length:int, type:str):
    if costIndex == 0:
        return 1
    elif costIndex == 1:
        return 1/coverage
    elif costIndex == 2: 
        return 1/math.sqrt(coverage)
    elif costIndex == 3:
        return i
        # if type == 'Exon':
        #     return i*(abs(length)-1)/coverage
        # else:
        #     return i/coverage

def flowDecompositionWithTranscriptlist(graph:dict, transcripts:list, decomposition_option:str, flow):
    #transcriptsCopy = deepcopy(transcripts)
    transcriptsCopy = list(nx.all_simple_paths(graph, '0', '1'))
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
        print('flowDecompositionDP with longestPath and residualFlow = ' + str(flow))

    elif decomposition_option == 'maximumFlow':
        # n = max(int(node) for node in list(graph.nodes()))+1
        # DP = np.zeros((n,n), dtype=int) # Initialize DP - Matrix
        # queue = [] # Initialize PriorityQueue
        # queue.append('0') # Add source to Priority Queue 
        # while(len(queue)>0): # Fill DP-Matrix
        #     v = queue.pop(0) 
        #     for u in list(graph.adj[v]):
        #         queue.append(u)
        #         DP[int(v)][int(u)] = graph.edges[v,u]['counts']['c']
        # flowSource = sum(DP[0][i] for i in range(len(DP)))
        # flowDrain = sum(DP[i][1] for i in range(len(DP)))
        # print(flowSource)
        # print(flowDrain)
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
            if max(DP[x][1] for x in range(len(DP)))==0:
                break 
            else:
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
                print(transcript)
        if flow>0:
            print(DP)
        print('flowDecompositionDP with maximumFlow and residualFlow = ' + str(flow))

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

