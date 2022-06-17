#!/usr/bin/env python

import sys, ast, os
from time import perf_counter
from xmlrpc.client import boolean

from sympy import true
import gurobipy as gp
import numpy as np
from collections import namedtuple 

ExonT = namedtuple('exon', 'exons count')
SpliceJunctionT = namedtuple('SpliceJunction', 'leftExon rightExon count')

def addCounts (graph:dict):
    exonSpliceJunctionCountDict = {}
    for edgeValue in graph.edges.values():                                                                  # Iterate over all edges in the Graph
        if edgeValue['type'] == 'Exon':
            exonSpliceJunctionCountDict[edgeValue['exon']] = edgeValue['counts']['c']
        elif edgeValue['type'] == 'SpliceJunction':
            exonSpliceJunctionCountDict[(edgeValue['startExon'], edgeValue['endExon'])] = edgeValue['counts']['c'] 
    return exonSpliceJunctionCountDict

def createAlphaIJMatrix (exonSpliceJunctionCounts:dict, localPathDict:dict):
    alpha=np.zeros((len(exonSpliceJunctionCounts), len(localPathDict)), dtype=int)
    i = 0
    for edge in exonSpliceJunctionCounts.keys():
        j = 0 
        for path in localPathDict.values():
            if type(edge) == tuple:
                for k in range(len(path)-1):
                    if path[k] == edge[0] and path[k+1] == edge[1]:
                        alpha[i][j]=1
            else:
                if edge in path:
                    alpha[i][j] = 1
            j = j+1
        i = i+1
    return alpha

def createAlphaIJDict(exonSpliceJunctionCounts:dict, localPathDict:dict):
    alpha={}
    for edge in exonSpliceJunctionCounts.keys():
        for path in localPathDict.values():
            alpha[(edge, tuple(path))]=0
            if type(edge) == tuple:
                for k in range(len(path)-1):
                    if path[k] == edge[0] and path[k+1] == edge[1]:
                        alpha[(edge, tuple(path))]=1
            else:
                if edge in path:
                    alpha[(edge, tuple(path))] = 1
    return alpha


def optimizeTranscripts (exonCounts:list, spliceJunctions:list, localPathDict:dict, alpha:np.array):
    m = gp.Model()
    
    # define data coefficients
    
    features = exonCounts + spliceJunctions
    transcripts = list(localPathDict.values())
    

    # add decision variables
    frequency = m.addVars(len(transcripts), vtype=gp.GRB.continuous, name='x')
    frequencySum = m.addVars(len(transcripts), vtype=gp.GRB.continuous, name='x')
    for j in range(len(features)):
        frequencySum[j] = gp.quicksum(frequency[i]*alpha[j][i] for i in range(len(frequency)))
    #
    formula = gp.quicksum(features[j].count - frequencySum[j] for j in range(len(features))) 
    m.setObjective(formula, gp.GRB.MINIMIZE)
    